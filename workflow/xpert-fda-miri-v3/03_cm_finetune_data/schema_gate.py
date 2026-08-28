#!/usr/bin/env python3
"""
Union-vs-v1 schema gate (step 03, primary seam).

Asserts that 03_cm_finetune_data/output/cm_union.h5ad is structure-identical
to v1's frozen finetune h5ad (workflow/xpert-fda-miri-v1/input/dtoxs.h5ad) in
every field the XPert library consumes, so the union trains through the
carried-over step scripts and the library without code changes:

  obs      — every v1 obs column used by the loader (pert_id, pert_idx,
             pert_dose, pert_dose_idx, pert_time, pert_time_idx, cell_iname,
             cell_idx, tissue_idx) present; dataset_id is an additive
             provenance column.
  var      — 978 landmarks in L1000 gene_idx order with gene_id (var_names,
             str Entrez) and gene_symbol.
  obsm     — "X_ctl" present, shape (n_obs, 978).
  X        — float32, log2(CPM+1) space (max in the v1-scale ~15 bound).

Any failure exits non-zero.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import (V1_DTOXS_H5AD, N_L1000_DRUGS,  # noqa: E402
                   AUTHOR_TIME_MAP, AUTHOR_TIMES, CM_DATASETS)

REQUIRED_OBS = ["pert_id", "pert_idx", "pert_dose", "pert_dose_idx",
                "pert_time", "pert_time_idx", "cell_iname", "cell_idx",
                "tissue_idx"]

# Union composition (ADR 0017): all four usable CM finetune datasets. The
# DToxS-family arms enter under v1's slot-24 convention (48 h in reality ->
# author 24 h slot, time_idx 2); gse244740/gse262419 are already 24 h. Every
# treated sample carries a positive real uM concentration (real-dose bins, ADR
# 0015 for the DToxS arms; the concentration-response datasets carry their own).
UNION_DATASETS = ["dtoxs_primary", "dtoxs_ipsc", "gse244740", "gse262419"]
QN_MAX = 15.002  # pretrain config max_value (l1000.yaml) / cm.yaml reset


def run(union_h5ad, v1_h5ad=None, verbose=True):
    import scanpy as sc

    v1_h5ad = v1_h5ad or V1_DTOXS_H5AD
    if not Path(v1_h5ad).exists():
        print(f"ERROR: v1 reference h5ad not found: {v1_h5ad}")
        sys.exit(1)

    a = sc.read_h5ad(union_h5ad)
    ref = sc.read_h5ad(v1_h5ad)
    errors = []

    # obs
    for col in REQUIRED_OBS:
        if col not in a.obs.columns:
            errors.append(f"obs column missing: {col}")
    if "dataset_id" not in a.obs.columns:
        errors.append("obs column missing: dataset_id (additive provenance)")

    # var — 978 landmarks in L1000 gene_idx order with gene_id / gene_symbol
    if a.n_vars != 978:
        errors.append(f"n_vars {a.n_vars} != 978")
    for col in ("gene_id", "gene_symbol"):
        if col not in a.var.columns:
            errors.append(f"var column missing: {col}")
    v1_var_ids = [str(v) for v in ref.var_names]
    union_var_ids = [str(v) for v in a.var_names]
    if v1_var_ids != union_var_ids:
        errors.append("var_names (gene_id order) differ from v1 (L1000 gene_idx order)")
    # var_names must equal gene_id strings (str Entrez)
    if list(a.var["gene_id"].astype(str)) != union_var_ids:
        errors.append("var.gene_id != var_names")

    # obsm
    if "X_ctl" not in a.obsm:
        errors.append("obsm['X_ctl'] missing")
    else:
        ctl = a.obsm["X_ctl"]
        if ctl.shape != (a.n_obs, 978):
            errors.append(f"obsm['X_ctl'] shape {ctl.shape} != ({a.n_obs}, 978)")

    # X dtype + log2 scale
    if str(a.X.dtype) != "float32":
        errors.append(f"X dtype {a.X.dtype} != float32")
    xmax = float(np.max(a.X.toarray())) if hasattr(a.X, "toarray") \
        else float(np.max(a.X))
    xmin = float(np.min(a.X.toarray())) if hasattr(a.X, "toarray") \
        else float(np.min(a.X))
    if xmax > 20.0:
        errors.append(f"X max {xmax:.2f} not in log2(CPM+1) scale (expected <= ~20)")

    # per-sample indices sane
    if a.obs["pert_idx"].min() < 0 or a.obs["pert_idx"].max() >= N_L1000_DRUGS:
        errors.append("pert_idx out of HG range [0, 8981)")
    if not a.obs["cell_idx"].isin([0]).all() or not a.obs["tissue_idx"].isin([0]).all():
        errors.append("cell_idx / tissue_idx must be 0 (all CM, num_tissue_id 1)")

    # --- author time/dose semantics (data-alignment gate) ---
    # Every sample must carry an author-valid time whose index equals
    # AUTHOR_TIME_MAP[time]. The DToxS-family datasets are encoded at the
    # author 24 h slot (v1 convention), so the union is all 24 h -> time_idx 2.
    bad_time = a.obs[~a.obs["pert_time"].isin(AUTHOR_TIMES)]
    if len(bad_time):
        errors.append(f"{len(bad_time)} samples with time outside {AUTHOR_TIMES}")
    bad_idx = a.obs[a.obs["pert_time"].map(AUTHOR_TIME_MAP) != a.obs["pert_time_idx"]]
    if len(bad_idx):
        errors.append(f"{len(bad_idx)} samples whose pert_time_idx != "
                      f"AUTHOR_TIME_MAP[pert_time]")
    # Zero dose is not allowed anywhere: every treated sample carries a real
    # positive uM concentration (ADR 0015; no single-dose fallback remains).
    bad_dose = a.obs[a.obs["pert_dose"] <= 0]
    if len(bad_dose):
        err_ds = sorted(bad_dose["dataset_id"].unique().tolist())
        errors.append(f"placeholder/zero dose in {err_ds} "
                      f"({len(bad_dose)} samples — every treated sample must "
                      f"carry a positive real concentration)")
    if not a.obs["dataset_id"].isin(UNION_DATASETS).all():
        err_ds = sorted(a.obs["dataset_id"].unique().tolist())
        errors.append(f"union composition {err_ds} not within {UNION_DATASETS}")

    # X within the pretrain bin range [0, QN_MAX] (config min/max reset)
    if xmax > QN_MAX:
        errors.append(f"X max {xmax:.2f} > QN_MAX {QN_MAX} (config bin range "
                      f"must be [0, {QN_MAX}])")
    if xmin < 0:
        errors.append(f"X min {xmin:.2f} < 0 (config bin range must be "
                      f"[0, {QN_MAX}])")

    # drug-inclusion gate: manifest n_union completeness + exclusions match
    # the union row counts per dataset (n_union counts only samples that
    # entered the union; n_samples/n_excl_* classify every extracted sample).
    manifest = pd.read_csv(Path(union_h5ad).parent / "drug_manifest.csv")
    if not manifest["hg_class"].isin(["in_hg", "extension", "unmatched"]).all():
        errors.append("drug manifest has unclassified compounds")
    for col in ("n_union", "n_excl_time", "n_excl_nodose"):
        if col not in manifest.columns:
            errors.append(f"drug manifest missing column: {col}")
    counts_by_ds = a.obs.groupby("dataset_id").size()
    man_in_hg = manifest[manifest["hg_class"] == "in_hg"] \
        .groupby("dataset")["n_union"].sum()
    for ds in counts_by_ds.index:
        if man_in_hg.get(ds, 0) != counts_by_ds[ds]:
            errors.append(f"dataset {ds}: union rows {counts_by_ds[ds]} != "
                          f"manifest in_hg n_union {man_in_hg.get(ds, 0)}")
    for ds in UNION_DATASETS:
        if ds not in counts_by_ds.index:
            errors.append(f"dataset {ds} absent from the union")
    dropped = set(CM_DATASETS) - set(counts_by_ds.index)
    for ds in dropped:
        if man_in_hg.get(ds, 0) != 0:
            errors.append(f"dataset {ds}: {man_in_hg.get(ds, 0)} in_hg n_union "
                          f"samples but no union rows (should be 0)")

    print(f"Schema gate — union {a.n_obs} x {a.n_vars} vs v1 {ref.n_obs} x {ref.n_vars}")
    for col in REQUIRED_OBS + ["dataset_id"]:
        print(f"  obs {col:15s} present")
    print(f"  var: 978 landmarks, gene_idx order matched: {v1_var_ids == union_var_ids}")
    print(f"  obsm X_ctl: {a.obsm.get('X_ctl').shape if 'X_ctl' in a.obsm else 'MISSING'}")
    print(f"  X dtype {a.X.dtype}, max {xmax:.2f}")
    print(f"  data-alignment: union {sorted(a.obs['dataset_id'].unique().tolist())}, "
          f"all time {sorted(a.obs['pert_time'].unique().tolist())} h "
          f"(time_idx {sorted(a.obs['pert_time_idx'].unique().tolist())}), "
          f"no zero-dose treated samples (real concentrations)")
    print(f"  drug-inclusion: pert_idx in [0, 8981) "
          f"({int(a.obs['pert_idx'].min())}-{int(a.obs['pert_idx'].max())}); "
          f"per-dataset union rows match manifest in_hg n_union counts")

    if errors:
        print("\nSCHEMA GATE FAILED:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    print("\nSchema gate passed: union is structure-identical to v1 in the "
          "fields the XPert loader consumes, with author time/dose semantics "
          "and pretrain-aligned expression values.")
    return True


if __name__ == "__main__":
    union = sys.argv[1] if len(sys.argv) > 1 \
        else str(Path(__file__).resolve().parent / "output" / "cm_union.h5ad")
    run(union)
