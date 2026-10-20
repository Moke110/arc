"""Bundle QA: cross-view invariants, Stage 1 regression and the QA report."""

from __future__ import annotations

import glob
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import paths
from .drugs import dataset_members
from .genes import GeneUniverse
from .sources import GSE244740_STAGE1_UNDERSCORE_FIXES, _dtoxs_key
from .views import ViewData

CONTRACT_OBS = [
    "well", "dataset_id", "compound_name", "drug_key", "inchikey", "smiles",
    "dose_um", "dose_source", "time_h", "is_control", "vehicle", "plate",
    "control_group", "cell_type", "donor", "replicate", "atc", "moa",
    "mesh_classes", "pharm_atc_moa", "pharm_mesh", "annotation_source",
    "lib_size", "n_genes_detected",
]


@dataclass
class GateResult:
    name: str
    passed: bool
    detail: str = ""
    skipped: bool = False


def view_exclusions(table, views: dict[str, ViewData]) -> pd.DataFrame:
    """Corpus compounds of a view's datasets that the view tier filters out."""
    rows = []
    for name, datasets, tier in paths.VIEWS:
        if name not in views:
            continue
        view = views[name]
        used = set(view.adata.obs.loc[~view.adata.obs["is_control"], "drug_key"])
        wanted = set(datasets)
        mask = table.compounds["datasets"].apply(
            lambda value: bool(dataset_members(value) & wanted))
        subset = table.compounds.loc[mask]
        for row in subset.to_dict("records"):
            if row["drug_key"] in used:
                continue
            rows.append({
                "view": name,
                "tier": tier,
                "drug_key": row["drug_key"],
                "compound_name": row["compound_name"],
                "datasets": row["datasets"],
                "reason": f"in corpus but outside the {tier} tier",
            })
    return pd.DataFrame(rows, columns=["view", "tier", "drug_key",
                                       "compound_name", "datasets", "reason"])


# ---------------------------------------------------------------------------
# Cross-view gates
# ---------------------------------------------------------------------------
def cross_view_gates(views: dict[str, ViewData],
                     universe: GeneUniverse) -> list[GateResult]:
    results = []

    # 1. identical var
    base = views[list(views)[0]].adata.var
    ok = all(view.adata.var.equals(base) for view in views.values())
    results.append(GateResult("var identical across all views", ok,
                              f"{len(base):,} genes"))

    # 2. gse262419 tier nesting + value identity + shared controls
    narrow_to_wide = [("gse262419_atc-moa.h5ad", "gse262419_atc-moa-mesh.h5ad"),
                      ("gse262419_atc-moa-mesh.h5ad", "gse262419_all-smiles.h5ad")]
    for narrow, wide in narrow_to_wide:
        if narrow not in views or wide not in views:
            continue
        a, b = views[narrow].adata, views[wide].adata
        a_pos = {well: i for i, well in enumerate(a.obs["well"])}
        b_pos = {well: i for i, well in enumerate(b.obs["well"])}
        missing = sorted(set(a_pos) - set(b_pos))
        ok = not missing
        detail = f"{narrow} rows absent from {wide}: {missing[:5]}"
        if ok:
            for well, row_a in a_pos.items():
                row_b = b_pos[well]
                if not np.array_equal(a.X[row_a], b.X[row_b]):
                    ok, detail = False, f"X differs for {well}"
                    break
                if not np.allclose(a.layers["log1p_cp10k"][row_a],
                                   b.layers["log1p_cp10k"][row_b], atol=1e-6):
                    ok, detail = False, f"layer differs for {well}"
                    break
                for column in CONTRACT_OBS:
                    if a.obs.iloc[row_a][column] != b.obs.iloc[row_b][column]:
                        ok, detail = False, f"{column} differs for {well}"
                        break
                if not ok:
                    break
        results.append(GateResult(f"tier nesting {narrow} ⊆ {wide}", ok, detail))

    # 3. the 272 gse262419 controls are identical across the three views
    gse_views = [views[name] for name in paths.GSE262419_VIEWS if name in views]
    if len(gse_views) >= 2:
        reference = gse_views[0].adata
        ref_ctrl = reference.obs["is_control"].to_numpy()
        ref_wells = reference.obs.loc[ref_ctrl, "well"].tolist()
        ok, detail = True, (f"{len(ref_wells)} control rows")
        for view in gse_views[1:]:
            adata = view.adata
            ctrl = adata.obs["is_control"].to_numpy()
            wells = adata.obs.loc[ctrl, "well"].tolist()
            if set(wells) != set(ref_wells):
                ok, detail = False, f"{view.name}: control wells differ"
                break
            position = {well: i for i, well in enumerate(adata.obs["well"])}
            for well in ref_wells:
                ref_row = reference.obs.index[reference.obs["well"] == well][0]
                ref_pos = reference.obs.index.get_loc(ref_row)
                row = position[well]
                if not np.array_equal(reference.X[ref_pos], adata.X[row]):
                    ok, detail = False, f"{view.name}: control counts differ at {well}"
                    break
                for column in ("control_group", "vehicle", "compound_name",
                               "dose_um"):
                    if (reference.obs.iloc[ref_pos][column]
                            != adata.obs.iloc[row][column]):
                        ok = False
                        detail = f"{view.name}: control {column} differs at {well}"
                        break
                if not ok:
                    break
            if not ok:
                break
        results.append(GateResult("gse262419 controls identical across views", ok, detail))
    return results


# ---------------------------------------------------------------------------
# Source-value spot checks (tickets 04/06/07)
# ---------------------------------------------------------------------------
def _raw_counts_value(dataset_id: str, well: str, gene_id: str) -> float | None:
    """Re-read one raw source cell for the spot check.

    Returns the value as the loader should see it (duplicate DToxS gene rows
    summed; duplicate gse262419 probe rows take the last occurrence, matching
    the loader).
    """
    raw_well = well.split("::", 1)[1]
    if dataset_id == "gse244740":
        frame = pd.read_csv(
            paths.GSE244740 / "GSE244740_processed_data_counts.txt.gz",
            sep="\t", usecols=["Ensembl", raw_well])
        values = frame.loc[frame["Ensembl"].astype(str) == gene_id, raw_well]
        return float(values.iloc[-1]) if len(values) else None

    if dataset_id == "gse262419":
        plate = raw_well.split("-")[0].replace("Plate", "")
        path = paths.GSE262419 / f"GSE262419_Plate{int(plate):02d}.csv.gz"
        frame = pd.read_csv(path)
        frame.columns = [str(column) for column in frame.columns]
        values = frame.loc[frame["Genes"].astype(str) == gene_id, raw_well]
        return float(values.iloc[-1]) if len(values) else None

    # DToxS arms: find the design row for the composite key, then the raw well.
    if dataset_id == "dtoxs_primary":
        designs = sorted(glob.glob(
            str(paths.DTOXS_LEVEL3 / "**"
                / "PromoCell-DGE-Experiment-Design-*.tsv"), recursive=True))
    elif dataset_id == "dtoxs_ipsc":
        designs = sorted(glob.glob(
            str(paths.DTOXS_IPSC / "GSE217423*Design*.tsv.gz")))
    else:
        return None
    for design_path in designs:
        design = pd.read_csv(design_path, sep="\t")
        keys = _dtoxs_key(design)
        match = keys == raw_well
        if not match.any():
            continue
        # dtoxs_primary columns are raw well ids; the ipsc DGE matrix already
        # uses the composite key as its column name.
        column = (design.loc[match, "Well"].iloc[0]
                  if dataset_id == "dtoxs_primary" else raw_well)
        counts_path = design_path.replace("Experiment-Design", "Read-Counts")
        try:
            frame = pd.read_csv(counts_path, sep="\t", index_col=0)
        except (ValueError, FileNotFoundError):
            continue
        values = frame.loc[frame.index.astype(str) == gene_id, column]
        return float(values.sum()) if len(values) else None
    return None


def source_value_check(sources: dict[str, object],
                       n_wells: int = 2) -> list[GateResult]:
    results = []
    for dataset_id, source in sources.items():
        meta = source.meta
        eligible = np.flatnonzero(
            (~meta["is_control"].to_numpy()) & (~meta["is_blank"].to_numpy()))
        checks, failures = [], []
        for row in eligible[:n_wells]:
            counts_row = source.counts[row]
            gene_position = int(np.argmax(counts_row))
            expected = int(counts_row[gene_position])
            if expected == 0:
                continue
            well = str(meta["well"].iloc[row])
            gene_id = str(source.gene_ids[gene_position])
            raw = _raw_counts_value(dataset_id, well, gene_id)
            checks.append((well, gene_id, raw, expected))
            if raw != expected:
                failures.append((well, gene_id, raw, expected))
        results.append(GateResult(
            f"source values {dataset_id}", not failures,
            (f"{len(checks)} sampled cells match the source file"
             if not failures else f"mismatches: {failures}"),
            skipped=not checks))
    return results


# ---------------------------------------------------------------------------
# Stage 1 regression (read-only)
# ---------------------------------------------------------------------------
def stage1_regression(sources: dict[str, object], universe: GeneUniverse,
                      n_wells: int = 5) -> list[GateResult]:
    import anndata as ad

    results = []
    for dataset_id in paths.SOURCE_DATASETS:
        path = paths.STAGE1_OUTPUT / f"{dataset_id}.h5ad"
        if not path.exists():
            results.append(GateResult(f"stage1 {dataset_id}", True,
                                      f"baseline missing: {path}", skipped=True))
            continue
        stage1 = ad.read_h5ad(path, backed="r")
        source = sources[dataset_id]
        stage1_wells = set(stage1.obs["well"].astype(str))
        our_wells = set(source.meta["well"].astype(str))
        missing = sorted(stage1_wells - our_wells)
        well_ok = not missing
        results.append(GateResult(
            f"stage1 {dataset_id}: well set", well_ok,
            (f"{len(stage1_wells):,} Stage 1 wells all present; "
             f"{len(our_wells - stage1_wells):,} extra source rows "
             "(controls / HG-excluded)") if well_ok else f"missing {missing[:5]}"))

        # Dose agreement on shared wells.  Stage 1 looked up underscore names
        # before normalization and fell back to the sunitinib top dose for the
        # four audited fixes; those wells are reported as corrections.
        shared = sorted(stage1_wells & our_wells)
        ours = source.meta.set_index("well")
        dose_mismatch, dose_fixed = [], []
        for well in shared:
            if well not in ours.index:
                continue
            expected = float(stage1.obs.loc[stage1.obs["well"] == well, "dose"].iloc[0])
            actual = float(ours.loc[well, "dose_um"])
            if abs(expected - actual) > 1e-9:
                name = str(ours.loc[well, "compound_name"])
                if name in GSE244740_STAGE1_UNDERSCORE_FIXES:
                    dose_fixed.append((well, expected, actual, name))
                else:
                    dose_mismatch.append((well, expected, actual, name))
        detail = f"{len(shared):,} shared wells"
        if dose_fixed:
            detail += (f"; {len(dose_fixed)} wells corrected from the Stage 1 "
                       f"underscore-name fallback "
                       f"({sorted({x[3] for x in dose_fixed})})")
        if dose_mismatch:
            detail = f"{len(dose_mismatch)} mismatches, e.g. {dose_mismatch[:2]}"
        results.append(GateResult(
            f"stage1 {dataset_id}: doses", not dose_mismatch, detail))

        # Counts: Stage 1 X is a per-plate monotonic (quantile) transform of
        # raw log2 CPM; ranks must agree exactly on the shared gene axis.
        # For gse262419, Stage 1's `entrez`-preferred landmark mapping selects
        # probes by the numeric suffix, which is a probe id, not an Entrez id;
        # it zeroes those landmarks.  Stage 2 keys on the symbol and recovers
        # them, so the rank check uses the Stage 1-mappable subset and reports
        # the exclusions explicitly.
        symbols = stage1.var["gene_symbol"].astype(str).tolist()
        allowed, excluded_entrez = set(symbols), 0
        stage1_map = (paths.STAGE1_OUTPUT / "gene_maps"
                      / f"{dataset_id}_to_l1000_gene_map.csv")
        if stage1_map.exists():
            mapping = pd.read_csv(stage1_map)
            if "match_type" in mapping:
                excluded_entrez = int((mapping["match_type"] == "entrez").sum())
                allowed = set(mapping.loc[
                    mapping["match_type"].isin(
                        ["direct", "legacy_alias", "ensembl"]),
                    "l1000_gene_symbol"].astype(str))
        universe_pos = {symbol: i for i, symbol in enumerate(universe.symbols)}
        affected = set()
        if dataset_id == "gse262419" and stage1_map.exists():
            # A probe suffix that equals *any* landmark Entrez is hijacked by
            # Stage 1's entrez-first lookup; Stage 2 keys on the symbol, so
            # those genes cannot be compared landmark-by-landmark.
            landmark_entrez = set(mapping["l1000_gene_id"].astype(str))
            symbol_of = dict(zip(universe.genes["hgnc_id"],
                                 universe.genes["symbol"]))
            for raw_id, hgnc_id in zip(source.gene_map.raw_ids,
                                       source.gene_map.hgnc_ids):
                if hgnc_id is None:
                    continue
                suffix = str(raw_id).rsplit("_", 1)[-1]
                if suffix in landmark_entrez:
                    affected.add(symbol_of.get(hgnc_id, ""))
        pairs = [(j, universe_pos[s]) for j, s in enumerate(symbols)
                 if s in universe_pos and s in allowed and s not in affected]
        if not pairs:
            results.append(GateResult(f"stage1 {dataset_id}: counts", False,
                                      "no shared genes"))
            continue
        columns, rows = zip(*pairs)
        sampled = shared[:n_wells]
        rhos = []
        for well in sampled:
            source_row = int(np.flatnonzero(source.meta["well"].to_numpy() == well)[0])
            raw = source.mapped_counts[source_row, list(rows)]
            stage1_row = int(np.flatnonzero(stage1.obs["well"].to_numpy() == well)[0])
            transformed = np.asarray(stage1.X[stage1_row]).ravel()[list(columns)]
            if np.all(raw == raw[0]):
                continue
            rho = float(spearmanr(raw, transformed).statistic)
            rhos.append(rho)
        ok = bool(rhos) and min(rhos) > 0.999
        detail = (f"{len(pairs)} shared genes, {len(rhos)} sampled wells, "
                  f"min Spearman {min(rhos):.6f}" if rhos
                  else f"{len(pairs)} shared genes, no testable wells")
        if excluded_entrez:
            detail += (f"; {excluded_entrez} Stage 1 entrez-mapped landmarks "
                       "excluded (probe-id suffix bug, recovered by symbol)")
        if affected:
            detail += (f"; {len(affected & allowed)} probe-suffix collision "
                       "genes excluded (Stage 1 hijacks the probe by suffix)")
        results.append(GateResult(f"stage1 {dataset_id}: counts rank", ok, detail))
    return results


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def _source_gene_stats(source, universe: GeneUniverse) -> dict:
    counts = dict(source.gene_map.counts)
    resolved = source.gene_map.resolved_ids
    return {
        "counts": counts,
        "resolved_unique": len(resolved),
        "entrez_verified": int(source.gene_map.verified.sum()),
        "in_intersection": len(resolved & set(universe.genes["hgnc_id"])),
    }


def write_report(path: Path, *, views: dict[str, ViewData],
                 sources: dict[str, object], universe: GeneUniverse,
                 table, gate_results: list[GateResult],
                 stage1_results: list[GateResult],
                 manifest: pd.DataFrame, drift_note: str) -> None:
    lines = [
        "# datasets-1 — QA report",
        "",
        f"Corpus: `corpus/{paths.CORPUS_VERSION}`",
        "",
        f"Gene axis: **{len(universe):,}** genes",
        "",
        f"Compounds in corpus: **{table.compounds['drug_key'].nunique():,}** drug keys",
        "",
        "## Per-view counts",
        "",
        "| view | perturbed | controls | drug keys | obs |",
        "|---|---|---|---|---|",
    ]
    for name, view in views.items():
        lines.append(f"| `{name}` | {view.n_perturbed:,} | {view.n_controls:,} | "
                     f"{view.n_drugs:,} | {view.adata.n_obs:,} |")
    lines += [
        "",
        "## Gates",
        "",
        "| gate | result | detail |",
        "|---|---|---|",
    ]
    for result in gate_results + stage1_results:
        status = ("SKIP" if result.skipped
                  else ("PASS" if result.passed else "FAIL"))
        lines.append(f"| {result.name} | {status} | {result.detail} |")

    lines += [
        "",
        "## Gene-side losses",
        "",
        "Per source: identifiers seen, how they resolved, and the unique genes",
        "they contributed. `unmapped`/`ambiguous` identifiers are dropped, never",
        "guessed.",
        "",
        "| source | direct | alias | previous | entrez | ensembl | entrez-verified | unmapped | ambiguous | unique genes | in intersection |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for dataset_id, source in sources.items():
        stats = _source_gene_stats(source, universe)
        counts = stats["counts"]
        lines.append(
            f"| {dataset_id} | {counts.get('direct', 0):,} | "
            f"{counts.get('alias', 0):,} | {counts.get('previous', 0):,} | "
            f"{counts.get('entrez', 0):,} | {counts.get('ensembl', 0):,} | "
            f"{stats['entrez_verified']:,} | {counts.get('unmapped', 0):,} | "
            f"{counts.get('ambiguous', 0):,} | {stats['resolved_unique']:,} | "
            f"{stats['in_intersection']:,} |")

    lines += [
        "",
        "## Excluded compounds (not in the corpus)",
        "",
    ]
    if table.exclusions.empty:
        lines.append("None.")
    else:
        lines += ["| dataset | compound | entity | wells | reason |",
                  "|---|---|---|---|---|"]
        for row in table.exclusions.to_dict("records"):
            lines.append(f"| {row['dataset_id']} | {row['compound_name']} | "
                         f"{row['entity_type']} | {row['n_wells']} | {row['reason']} |")

    lines += [
        "",
        "## View-level exclusions (in the corpus, filtered by the view tier)",
        "",
    ]
    excluded = view_exclusions(table, views)
    for name, datasets, tier in paths.VIEWS:
        if name not in views:
            continue
        group = excluded[excluded["view"] == name] if not excluded.empty else excluded
        if group.empty:
            lines.append(f"- `{name}` ({tier}): no view-level compound exclusions")
            continue
        entries = len(group)
        keys = group["drug_key"].nunique()
        lines.append(f"<details><summary><code>{name}</code> ({tier}): "
                     f"{entries} compound {'entry' if entries == 1 else 'entries'} "
                     f"({keys} drug {'key' if keys == 1 else 'keys'}) filtered "
                     f"by the tier</summary>")
        lines.append("")
        lines.append("| compound | drug_key | reason |")
        lines.append("|---|---|---|")
        for row in group.sort_values("compound_name").to_dict("records"):
            lines.append(f"| {row['compound_name']} | `{row['drug_key']}` | "
                         f"{row['reason']} |")
        lines.append("")
        lines.append("</details>")
        lines.append("")

    lines += [
        "",
        "## Salt / free-base collapses",
        "",
        "These names share one `drug_key` (parent structure) but remain separate",
        "rows, so a cold-drug split cannot leak across them:",
        "",
    ]
    pairs = [("Tamoxifen", "Tamoxifen citrate"),
             ("Diltiazem", "Diltiazem hydrochloride"),
             ("Estradiol", "17beta-Estradiol"),
             ("Perfluorooctanoic acid", "Ammonium perfluorooctanoate")]
    by_name = table.compounds.set_index("compound_name")
    for free, salt in pairs:
        if free in by_name.index and salt in by_name.index:
            lines.append(f"- {free} / {salt}: `{by_name.loc[free, 'drug_key']}`")
    lines += [
        "",
        "Doxorubicin and epirubicin are 4'-epimers and carry **different** parent",
        "keys (connectivity-only keys would collapse them; the build rejects",
        "14-character keys).",
        "",
        "## Tier counts",
        "",
        "| tier | drug keys (corpus) |",
        "|---|---|",
    ]
    for tier in ("atc-moa", "atc-moa-mesh", "all-smiles"):
        lines.append(f"| {tier} | {len(table.members(tier)):,} |")
    lines += [
        "",
        "## Drift",
        "",
        drift_note,
        "",
        "## Manifest",
        "",
        "| view | obs | perturbed | controls | drug keys | genes | content sha256 |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in manifest.to_dict("records"):
        lines.append(f"| `{row['view']}` | {row['n_obs']:,} | "
                     f"{row['n_perturbed']:,} | {row['n_controls']:,} | "
                     f"{row['n_drugs']:,} | {row['n_genes']:,} | "
                     f"`{row['content_sha256'][:16]}` |")
    lines += [
        "",
        "## Caveats",
        "",
        "- `dtoxs_primary` doses are literature/clinical Cmax values; `DIC`",
        "  (diclofenac, 3 wells) uses the 5 uM recorded in the DToxS DCIC design",
        "  workbook and is marked `dose_source=measured`.",
        "- `gse244740` doses are derived from the curated paper Table-1 top",
        "  concentration (`dose_source=derived_top_dose`), not measured.",
        "- The DToxS `AMI` arm is a pentyl-benzofuran homolog (CID 16218975), not",
        "  amiodarone; it is in `all-smiles` only.",
        "- `dtoxs_cm.h5ad` intentionally mixes `primary_cm` and `ipsc_cm`; cell",
        "  type and donor are covariates, not collapsed constants.",
        "- The 32 `gse262419` `Blank` wells are technical blanks, neither",
        "  controls nor perturbations.",
        "",
    ]
    path.write_text("\n".join(lines))
