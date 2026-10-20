"""Freeze and verify ``corpus/v1/`` (gene set + drug set, no matrices).

The corpus is the only thing that makes a rerun comparable: a bundle re-reads
the sources but applies the frozen universes.  Any change to the computed
intersection or to a resolved drug key is a **hard failure**; promoting a new
corpus version is a deliberate edit (ADR 0026 decision 6).
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import pandas as pd

from . import paths
from .drugs import CompoundTable
from . import genes
from .genes import GeneUniverse


class DriftError(RuntimeError):
    """A rerun produced a gene or drug outside the frozen corpus."""


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")


def freeze_corpus(universe: GeneUniverse, table: CompoundTable,
                  source_stats: dict, source_maps: dict, index,
                  out_dir: Path = paths.CORPUS_DIR,
                  rebuild: bool = False) -> None:
    if out_dir.exists() and not rebuild:
        raise FileExistsError(f"corpus already exists: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    genes = universe.genes.copy()
    _write_csv(genes, out_dir / "genes.csv")
    _write_csv(table.compounds, out_dir / "compounds.csv")
    _write_csv(table.exclusions, out_dir / "excluded_compounds.csv")

    # One row per (source, resolved HGNC gene).  This is the frozen per-source
    # gene set: a rerun that maps a new gene in any source fails, even when the
    # gene is outside the four-way intersection and would otherwise be silently
    # dropped by the view writer.
    source_rows = []
    for dataset_id, mapping in source_maps.items():
        seen = set()
        for hgnc_id, match, verified in zip(mapping.hgnc_ids, mapping.match_types,
                                            mapping.verified):
            if hgnc_id is None or hgnc_id in seen:
                continue
            seen.add(hgnc_id)
            source_rows.append({
                "dataset_id": dataset_id,
                "hgnc_id": hgnc_id,
                "symbol": index.ref(hgnc_id).symbol,
                "match_type": match,
                "entrez_verified": bool(verified),
            })
    _write_csv(pd.DataFrame(source_rows).sort_values(
        ["dataset_id", "symbol"], kind="stable"), out_dir / "source_gene_sets.csv")

    tier_counts = {
        tier: table.members(tier)
        for tier in ("atc-moa", "atc-moa-mesh", "all-smiles")
    }
    match_counts = genes["match_type"].value_counts().to_dict()
    lines = [
        "# corpus/v1 — frozen Stage 2 CM gene and drug universes",
        "",
        f"**Frozen:** {date.today().isoformat()}",
        "",
        f"**Source commit (vendored Stage 1 loaders):** `{paths.SOURCE_COMMIT}`",
        "",
        f"**Genes:** {len(genes):,}",
        "",
        f"**Compounds:** {table.compounds['drug_key'].nunique():,} drug keys "
        f"({len(table.compounds):,} key/name rows)",
        "",
        "## Gene axis",
        "",
        "Four-source intersection (`dtoxs_primary`, `dtoxs_ipsc`, `gse244740`,",
        "`gse262419`). Each source identifier resolves to exactly one current",
        "HGNC symbol:",
        "",
        "1. DToxS symbols are Entrez-verified through the package's own",
        "   `UCSC-hg38-Gene-HGNC-Entrez-Map.tsv` (an Entrez ID wins when HGNC",
        "   has reassigned the symbol).",
        "2. `gse244740` Ensembl IDs (version stripped) map through",
        "   `ensembl_gene_id`.",
        "3. `gse262419` `SYMBOL_probeid` rows are keyed by symbol; duplicate",
        "   probes collapse to one gene and the numeric suffix is only an",
        "   Entrez fallback (it is a probe id, not an Entrez id).",
        "",
        "Unresolvable or ambiguous identifiers are dropped and counted in",
        "`../README.md` and the bundle `qa_report.md`; nothing is guessed.",
        "The order is the HGNC symbol, lexicographic, and is authoritative for",
        "every bundle `var`.",
        "",
        f"Match types: {json.dumps(match_counts, sort_keys=True)}",
        "",
        "## Drug set",
        "",
        "`drug_key` is the desalted parent InChIKey **with stereochemistry**.",
        "`inchikey` is the exact registered structure key, so salt and free-base",
        "rows share a `drug_key` but keep distinct `inchikey` values. The three",
        "tiers are boolean columns:",
        "",
        "- `pharm_atc_moa`: ATC code or ChEMBL mechanism;",
        "- `pharm_mesh`: the above or a PubChem MeSH class;",
        "- `pharm_all_smiles`: a resolvable small-molecule structure.",
        "",
        f"Tier sizes: atc-moa {len(tier_counts['atc-moa'])}, "
        f"atc-moa-mesh {len(tier_counts['atc-moa-mesh'])}, "
        f"all-smiles {len(tier_counts['all-smiles'])} drug keys.",
        "",
        "## Inputs",
        "",
        "- HGNC complete set: `pert_data/annotation/hgnc_complete_set.txt`",
        "- compound registry snapshot:",
        "  `pert_data/annotation/cm_compound_annotation.csv`",
        "- source data: `pert_data/{dtoxs_primary,dtoxs_ipsc,gse244740,gse262419}`",
        "",
        "Both annotation assets are pinned with URL, retrieval date and",
        "checksum in `pert_data/annotation/README.md`; the build makes no",
        "network call.",
        "",
        "`source_gene_sets.csv` is the frozen per-source resolved gene set: a",
        "rerun that maps a new gene in any source fails, even if the gene is",
        "outside the four-way intersection.",
        "",
        "## Per-source gene resolution",
        "",
        "| source | direct | alias | previous | entrez | ensembl | unmapped | ambiguous | resolved |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for dataset, stats in source_stats.items():
        counts = stats["counts"]
        resolved = stats["resolved_unique"]
        lines.append(
            f"| {dataset} | {counts.get('direct', 0):,} | "
            f"{counts.get('alias', 0):,} | {counts.get('previous', 0):,} | "
            f"{counts.get('entrez', 0):,} | {counts.get('ensembl', 0):,} | "
            f"{counts.get('unmapped', 0):,} | {counts.get('ambiguous', 0):,} | "
            f"{resolved:,} |")
    lines += [
        "",
        "See the bundle `gene_crosswalk.csv` for the raw identifier of every",
        "gene and the bundle `qa_report.md` for per-view gates.",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines))


def load_frozen_corpus(out_dir: Path = paths.CORPUS_DIR
                       ) -> tuple[GeneUniverse, pd.DataFrame, pd.DataFrame]:
    if not (out_dir / "genes.csv").exists():
        raise FileNotFoundError(f"frozen corpus missing: {out_dir / 'genes.csv'}")
    universe = genes.universe_from_genes_csv(out_dir / "genes.csv")
    compounds = pd.read_csv(out_dir / "compounds.csv")
    source_sets = pd.read_csv(out_dir / "source_gene_sets.csv")
    return universe, compounds, source_sets


def check_drift(computed: GeneUniverse, table: CompoundTable,
                source_maps: dict, index,
                frozen: GeneUniverse, frozen_compounds: pd.DataFrame,
                frozen_source_sets: pd.DataFrame) -> None:
    added = sorted(set(computed.symbols) - set(frozen.symbols))
    removed = sorted(set(frozen.symbols) - set(computed.symbols))
    source_issues = []
    for dataset_id, mapping in source_maps.items():
        computed_ids = {h for h in mapping.hgnc_ids if h is not None}
        frozen_ids = set(frozen_source_sets.loc[
            frozen_source_sets["dataset_id"] == dataset_id, "hgnc_id"])
        new = sorted(computed_ids - frozen_ids)
        lost = sorted(frozen_ids - computed_ids)
        if new or lost:
            new_symbols = [index.ref(h).symbol for h in new[:10]]
            lost_symbols = frozen_source_sets.loc[
                frozen_source_sets["hgnc_id"].isin(lost), "symbol"].tolist()[:10]
            source_issues.append(
                f"  {dataset_id}: added {len(new)} {new_symbols}, "
                f"removed {len(lost)} {lost_symbols}")
    new_pairs = set(zip(table.compounds["drug_key"],
                        table.compounds["compound_name"]))
    frozen_pairs = set(zip(frozen_compounds["drug_key"],
                           frozen_compounds["compound_name"]))
    new_keys = sorted({key for key, _ in new_pairs - frozen_pairs})
    lost_keys = sorted({key for key, _ in frozen_pairs - new_pairs})
    if added or removed or new_keys or lost_keys or source_issues:
        lines = ["corpus drift detected (a deliberate corpus version bump is required):"]
        if added:
            lines.append(f"  genes added ({len(added)}): {added[:10]}")
        if removed:
            lines.append(f"  genes removed ({len(removed)}): {removed[:10]}")
        if new_keys:
            lines.append(f"  drug key/name pairs added ({len(new_keys)}): {new_keys[:10]}")
        if lost_keys:
            lines.append(f"  drug key/name pairs removed ({len(lost_keys)}): {lost_keys[:10]}")
        lines.extend(source_issues)
        raise DriftError("\n".join(lines))
