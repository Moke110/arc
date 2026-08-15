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
from utils import V1_DTOXS_H5AD, N_L1000_DRUGS  # noqa: E402

REQUIRED_OBS = ["pert_id", "pert_idx", "pert_dose", "pert_dose_idx",
                "pert_time", "pert_time_idx", "cell_iname", "cell_idx",
                "tissue_idx"]


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
    if xmax > 20.0:
        errors.append(f"X max {xmax:.2f} not in log2(CPM+1) scale (expected <= ~20)")

    # per-sample indices sane
    if a.obs["pert_idx"].min() < 0 or a.obs["pert_idx"].max() >= N_L1000_DRUGS:
        errors.append("pert_idx out of HG range [0, 8981)")
    if not a.obs["pert_time_idx"].isin([0, 1, 2]).all():
        errors.append("pert_time_idx not in union map {24:0, 48:1, 144:2}")
    if not a.obs["cell_idx"].isin([0]).all() or not a.obs["tissue_idx"].isin([0]).all():
        errors.append("cell_idx / tissue_idx must be 0 (all CM, num_tissue_id 1)")

    # drug-inclusion gate: manifest completeness + exclusions match row counts
    manifest = pd.read_csv(Path(union_h5ad).parent / "drug_manifest.csv")
    if not manifest["hg_class"].isin(["in_hg", "extension", "unmatched"]).all():
        errors.append("drug manifest has unclassified compounds")
    counts_by_ds = a.obs.groupby("dataset_id").size()
    man_in_hg = manifest[manifest["hg_class"] == "in_hg"] \
        .groupby("dataset")["n_samples"].sum()
    for ds in counts_by_ds.index:
        if man_in_hg.get(ds, 0) != counts_by_ds[ds]:
            errors.append(f"dataset {ds}: union rows {counts_by_ds[ds]} != "
                          f"manifest in_hg samples {man_in_hg.get(ds, 0)}")

    print(f"Schema gate — union {a.n_obs} x {a.n_vars} vs v1 {ref.n_obs} x {ref.n_vars}")
    for col in REQUIRED_OBS + ["dataset_id"]:
        print(f"  obs {col:15s} present")
    print(f"  var: 978 landmarks, gene_idx order matched: {v1_var_ids == union_var_ids}")
    print(f"  obsm X_ctl: {a.obsm.get('X_ctl').shape if 'X_ctl' in a.obsm else 'MISSING'}")
    print(f"  X dtype {a.X.dtype}, max {xmax:.2f}")
    print(f"  drug-inclusion: pert_idx in [0, 8981) "
          f"({int(a.obs['pert_idx'].min())}-{int(a.obs['pert_idx'].max())}); "
          f"per-dataset union rows match manifest in_hg counts")

    if errors:
        print("\nSCHEMA GATE FAILED:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    print("\nSchema gate passed: union is structure-identical to v1 in the "
          "fields the XPert loader consumes.")
    return True


if __name__ == "__main__":
    union = sys.argv[1] if len(sys.argv) > 1 \
        else str(Path(__file__).resolve().parent / "output" / "cm_union.h5ad")
    run(union)
