"""Corpus QA gates for the native-axis Stage 2 CM corpora (ticket #15).

Two independent checks:

* ``corpus_structure_gates`` runs **in the build** against the in-memory
  ``AnnData`` before it is written: gate G1 (the h5ad contract) plus the row /
  chemical counts of the published corpus table (G2 counting).
* ``source_value_check`` re-reads a raw source cell for a couple of wells and
  compares it with the loaded matrix, so a silent loader change cannot reach a
  corpus unnoticed.

The retired v1 view gates (shared var axis, tier nesting, per-view
``control_mean``) were removed with the view layer (ADR 0033 decision 5); the
frozen copies live in ``stage2/_frozen/cm-data-v1/execution/``.
"""

from __future__ import annotations

import glob
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import paths
from .sources import _dtoxs_key


@dataclass
class GateResult:
    name: str
    passed: bool
    detail: str = ""
    skipped: bool = False


def _fail(results: list[GateResult], name: str, detail: str) -> None:
    results.append(GateResult(name, False, detail))


def corpus_structure_gates(adata, corpus_id: str,
                           expected: dict | None = None) -> list[GateResult]:
    """G1 structure contract + G2 counts for one built corpus."""
    results: list[GateResult] = []
    obs, var, X = adata.obs, adata.var, adata.X

    # 1. unique row index (one row per raw well / profile)
    wells = obs.index.astype(str)
    ok = bool(obs.index.is_unique and not obs.index.isna().any()
              and (wells != "").all())
    results.append(GateResult(
        "index unique", ok,
        f"{len(wells):,} rows" if ok else "duplicate, empty or missing row index"))

    # 2. obs contract columns, in order, no integer encodings
    missing = [column for column in paths.CORPUS_OBS
               if column not in obs.columns]
    if missing:
        _fail(results, "obs contract columns", f"missing {missing}")
    else:
        results.append(GateResult(
            "obs contract columns", True,
            f"{len(paths.CORPUS_OBS)} columns in the ADR-0033 order"))
    forbidden = {"pert_idx", "cell_idx", "dose_idx", "pert_dose_idx"}
    leaked = sorted(forbidden & set(obs.columns))
    results.append(GateResult(
        "no integer model encodings", not leaked,
        "none" if not leaked else f"leaked columns {leaked}"))

    # 3. X = non-negative integer raw counts
    integer = np.issubdtype(X.dtype, np.integer)
    finite = bool(np.isfinite(X).all()) if np.issubdtype(
        X.dtype, np.floating) else True
    negative = int((X < 0).sum())
    ok = bool(integer and finite and negative == 0)
    results.append(GateResult(
        "X non-negative integer counts", ok,
        (f"{X.dtype}, shape {X.shape}, min 0" if ok
         else f"dtype {X.dtype}, {negative} negative cells, finite={finite}")))

    # 4. the log1p_cp10k layer equals log1p(1e4 * x / row_sum) (chunked)
    lib_size = X.sum(axis=1, dtype=np.int64)
    layer_ok, layer_detail = False, "layer missing"
    layer = adata.layers.get("log1p_cp10k")
    if layer is not None:
        layer = np.asarray(layer)
        layer_ok = True
        for start in range(0, X.shape[0], 512):
            stop = min(start + 512, X.shape[0])
            expected_layer = np.log1p(X[start:stop] / lib_size[start:stop, None]
                                      * 1e4)
            if not np.allclose(layer[start:stop], expected_layer,
                               rtol=1e-5, atol=1e-6):
                layer_ok, layer_detail = False, f"differs in rows {start}..{stop}"
                break
        if layer_ok:
            layer_detail = (f"float32 formula holds on {X.shape[0]:,} rows"
                            if layer.dtype == np.float32
                            else f"formula holds but dtype is {layer.dtype}")
            layer_ok = bool(layer.dtype == np.float32)
    results.append(
        GateResult("log1p_cp10k layer", layer_ok, layer_detail))

    # 5. no extra containers / uns payload
    uns = {"has": list(getattr(adata, "uns", {}).keys())}
    uns_ok = not uns["has"]
    results.append(GateResult(
        "uns empty", uns_ok, "empty" if uns_ok else f"keys {uns['has']}"))

    containers = {name: list(getattr(adata, name, {}).keys())
                  for name in ("obsm", "varm", "obsp", "varp")}
    extra = {name: keys for name, keys in containers.items() if keys}
    results.append(GateResult(
        "obsm/varm/obsp/varp empty", not extra,
        "empty" if not extra else f"present {sorted(extra)}"))

    # 6. per-row library statistics
    n_detected = (X > 0).sum(axis=1).astype(np.int32)
    ok = bool((obs["lib_size"].to_numpy() == lib_size).all()
              and (obs["n_genes_detected"].to_numpy() == n_detected).all()
              and (lib_size > 0).all())
    results.append(GateResult(
        "lib_size / n_genes_detected", ok,
        f"row sums match on {X.shape[0]:,} rows"
        if ok else "obs library statistics do not match X"))

    # 7. row semantics: controls, blanks, perturbed compounds
    control = obs["is_control"].to_numpy()
    blank = obs["is_blank"].to_numpy()
    perturbed = ~control & ~blank
    problems = []
    if (obs.loc[control, "dose_um"].to_numpy() != 0).any():
        problems.append("control rows with a non-zero dose")
    if (obs.loc[control, "drug_key"].astype(str) != "").any():
        problems.append("control rows with a drug_key")
    if (blank & control).any():
        problems.append("rows both blank and control")
    if (obs.loc[blank, "dose_um"].to_numpy() != 0).any():
        problems.append("blank rows with a non-zero dose")
    if (obs.loc[perturbed, "dose_um"].to_numpy() <= 0).any():
        problems.append("perturbed rows without a positive dose")
    for column in ("drug_key", "smiles", "inchikey"):
        if (obs.loc[perturbed, column].astype(str) == "").any():
            problems.append(f"perturbed rows missing {column}")
    results.append(GateResult(
        "row semantics", not problems,
        "controls/blanks/perturbed rows well formed" if not problems
        else "; ".join(problems)))

    # 8. single dataset identity + published counts
    dataset_ids = set(obs["dataset_id"].astype(str))
    ok = dataset_ids == {corpus_id}
    stats = {
        "rows": int(adata.n_obs),
        "perturbed": int(perturbed.sum()),
        "control": int(control.sum()),
        "blank": int(blank.sum()),
        "chemicals": int(obs.loc[perturbed, "drug_key"].astype(str).nunique()),
        "axis": int(adata.n_vars),
    }
    results.append(GateResult(
        "dataset_id", ok,
        f"{{'{corpus_id}'}}" if ok else f"got {sorted(dataset_ids)}"))
    if expected:
        mismatched = {key: (expected[key], stats.get(key))
                      for key in expected if stats.get(key) != expected[key]}
        results.append(GateResult(
            "published counts", not mismatched,
            f"rows {stats['rows']:,} / pert {stats['perturbed']:,} / "
            f"ctrl {stats['control']:,} / blank {stats['blank']:,} / "
            f"chemicals {stats['chemicals']:,} / axis {stats['axis']:,}"
            if not mismatched else f"mismatched {mismatched}"))
    return results


# ---------------------------------------------------------------------------
# Source-value spot checks
# ---------------------------------------------------------------------------
def _raw_counts_value(dataset_id: str, well: str, gene_id: str) -> float | None:
    """Re-read one raw source cell for the spot check.

    Returns the value as the loader sees it (duplicate DToxS gene rows summed;
    duplicate temp_env probe rows take the last occurrence, matching the
    loader).
    """
    raw_well = well.split("::", 1)[1]
    if dataset_id == "screen_ctox":
        frame = pd.read_csv(
            paths.SCREEN_CTOX_DIR / "GSE244740_processed_data_counts.txt.gz",
            sep="\t", usecols=["Ensembl", raw_well])
        values = frame.loc[frame["Ensembl"].astype(str) == gene_id, raw_well]
        return float(values.iloc[-1]) if len(values) else None

    if dataset_id == "temp_env":
        plate = raw_well.split("-")[0].replace("Plate", "")
        path = paths.TEMP_ENV_DIR / f"GSE262419_Plate{int(plate):02d}.csv.gz"
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
