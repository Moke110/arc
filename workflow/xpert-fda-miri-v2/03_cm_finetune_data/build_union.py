#!/usr/bin/env python3
"""
Step 03 — CM finetune-data union (the substantive v2 change).

Reads all five CM datasets under pert_data/cm **read-only** and assembles the
workflow-owned union finetune h5ad that steps 05/06/07 train on:

  03_cm_finetune_data/output/cm_union.h5ad

The union h5ad is schema-identical to v1's frozen input/dtoxs.h5ad in the
fields the XPert library consumes (obs pert_id/pert_idx/pert_dose[_idx]/
pert_time[_idx]/cell_iname/cell_idx/tissue_idx + the additive dataset_id
provenance column; var = the 978 L1000 landmarks in gene_idx order with
gene_id/gene_symbol; obsm["X_ctl"] per-sample matched control; X float32 in
log2(CPM+1) space). Only samples whose compound resolves to an HG drug
(drug_idx in [0, 8981)) enter the union.

Also emits:
  03_cm_finetune_data/output/gene_maps/<dataset>_to_l1000_gene_map.csv
      per-dataset landmark resolution (direct | legacy_alias | entrez | ensembl)
  03_cm_finetune_data/output/drug_manifest.csv
      per-compound drug resolution (in_hg | extension | unmatched)
  03_cm_finetune_data/output/dataset_summary.csv
      per-dataset sample / resolution / control counts

Usage:
    python 03_cm_finetune_data/build_union.py            # idempotent
    python 03_cm_finetune_data/build_union.py --force    # rebuild
    python 03_cm_finetune_data/build_union.py --skip-schema-gate
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import _datasets  # noqa: E402
import gene_mapping  # noqa: E402
from drug_resolution import build_drug_manifest  # noqa: E402
from utils import (WORKFLOW_ROOT, CM_DATASETS, UNION_TIME_MAP,  # noqa: E402
                   N_L1000_DRUGS)

STEP_DIR = WORKFLOW_ROOT / "03_cm_finetune_data"
OUT_DIR = STEP_DIR / "output"
GENE_MAP_DIR = OUT_DIR / "gene_maps"
UNION_H5AD = OUT_DIR / "cm_union.h5ad"
MANIFEST = OUT_DIR / "drug_manifest.csv"
SUMMARY = OUT_DIR / "dataset_summary.csv"
CELL_NAME = "cardio"


def dose_idx_of(dose):
    """XPert assign_pert_dose_corrected bins (0-9) — mirrored from the library."""
    if 0 <= dose < 0.21:
        return 0
    if dose < 0.41:
        return 1
    if dose < 0.71:
        return 2
    if dose < 1.01:
        return 3
    if dose < 1.51:
        return 4
    if dose < 3.1:
        return 5
    if dose < 4.1:
        return 6
    if dose < 7.1:
        return 7
    if dose < 12.1:
        return 8
    return 9


def build_gene_maps():
    maps = {}
    for ds in CM_DATASETS:
        mp = gene_mapping.build_dataset_map(ds)
        gene_mapping.write_dataset_map(ds, mp, GENE_MAP_DIR)
        maps[ds] = mp
        n = (mp["match_type"] != "unmatched").sum()
        print(f"  {ds}: {n}/978 landmarks resolved "
              f"({dict(mp['match_type'].value_counts())})")
    return maps


def source_to_landmark_col(gene_map):
    """Row-source identifier -> landmark column index (0..977).

    For 'symbol' / 'ensembl' spaces the source is the row identifier itself;
    for 'symbol_entrez' the source is either the numeric part (entrez matches)
    or the symbol part (direct / legacy_alias matches).
    """
    meta = gene_mapping.load_l1000_meta()
    idx_by_gid = {str(g): i for i, g in zip(meta["gene_idx"], meta["gene_id"])}
    out = {}
    for _, r in gene_map.iterrows():
        src = r["source_symbol"]
        if src:
            out[src] = idx_by_gid[str(r["l1000_gene_id"])]
    return out


def build_landmark_matrix(counts, gene_ids, gene_map, space):
    """log2(CPM+1) (n_samples, 978) matrix in L1000 gene_idx order.

    Duplicate source rows are summed (v1 convention); source identifiers that
    do not map to a landmark are dropped.
    """
    meta = gene_mapping.load_l1000_meta()
    idx_by_gid = {str(g): i for i, g in zip(meta["gene_idx"], meta["gene_id"])}
    row2col = source_to_landmark_col(gene_map)
    if space == "symbol_entrez":
        mt_by_src = {r["source_symbol"]: r["match_type"]
                     for _, r in gene_map.iterrows() if r["source_symbol"]}
        row2col_ent = {src: c for src, c in row2col.items()
                       if mt_by_src.get(src) == "entrez"}
        row2col_sym = {src: c for src, c in row2col.items()
                       if mt_by_src.get(src) != "entrez"}
    # map each gene row
    row_col = np.full(len(gene_ids), -1, dtype=int)
    for i, gid in enumerate(gene_ids):
        if space == "symbol_entrez":
            parts = str(gid).rsplit("_", 1)
            num = parts[1] if len(parts) == 2 else ""
            if num in row2col_ent:
                row_col[i] = row2col_ent[num]
            elif parts[0] in row2col_sym:
                row_col[i] = row2col_sym[parts[0]]
        else:
            key = str(gid).split(".")[0] if space == "ensembl" else str(gid)
            row_col[i] = row2col.get(key, -1)

    # sum duplicate rows per landmark column
    mat = np.zeros((counts.shape[0], 978), dtype=np.float64)
    valid = row_col >= 0
    np.add.at(mat.T, row_col[valid], counts[:, valid].astype(np.float64).T)
    # CPM + log2(CPM+1) per sample
    totals = mat.sum(axis=1)
    cpm = np.divide(mat, totals[:, None], out=np.zeros_like(mat),
                    where=totals[:, None] > 0) * 1e6
    return np.log2(cpm + 1).astype(np.float32)


def matched_control_matrix(meta, X_landmark, controls_mask, group_col):
    """Per-sample matched control = group mean of control log2(CPM+1)."""
    ctl = np.zeros_like(X_landmark, dtype=np.float64)
    if controls_mask.sum() == 0:
        return ctl.astype(np.float32)
    groups = meta.loc[controls_mask, group_col].unique()
    for g in groups:
        sel = meta[group_col] == g
        if controls_mask[sel.values].sum() == 0:
            continue
        ctl_rows = X_landmark[controls_mask & sel.values]
        grp_mean = ctl_rows.mean(axis=0)
        for i in np.where(~controls_mask & sel.values)[0]:
            ctl[i] = grp_mean
    return ctl.astype(np.float32)


def process_dataset(ds, gene_map, drug_idx_of):
    meta, counts, gene_ids = _datasets.extract_samples(ds)
    space = gene_mapping.ID_SPACE[ds]
    X = build_landmark_matrix(counts, gene_ids, gene_map, space)

    # matched controls
    if ds in ("gse262419", "gse244740"):
        grp_col = "plate"
    else:
        grp_col = "dataset_id"
    X_ctl = matched_control_matrix(meta, X, meta["is_control"].values, grp_col)

    # restrict to HG-resolved treated samples
    trt = ~meta["is_control"].values
    drug_idx = np.array([drug_idx_of.get((ds, str(r.compound_code) if r.compound_code else r.compound_name), -1)
                         for r in meta.itertuples()])
    keep = trt & (drug_idx >= 0) & (drug_idx < N_L1000_DRUGS)

    obs = meta[keep].reset_index(drop=True).copy()
    obs["pert_idx"] = drug_idx[keep]
    obs["pert_dose"] = obs["dose"]
    obs["pert_dose_idx"] = obs["dose"].map(dose_idx_of)
    obs["pert_time"] = obs["time"]
    obs["pert_time_idx"] = obs["time"].map(UNION_TIME_MAP)
    obs["cell_iname"] = CELL_NAME
    obs["cell_idx"] = 0
    obs["tissue_idx"] = 0
    obs["pert_id"] = obs["well"]
    return obs, X[keep], X_ctl[keep]


def main():
    parser = argparse.ArgumentParser(description="Build the CM finetune-data union")
    parser.add_argument("--force", action="store_true",
                        help="Rebuild even if the union h5ad exists")
    parser.add_argument("--skip-schema-gate", action="store_true",
                        help="Skip the union-vs-v1 schema comparison")
    args = parser.parse_args()

    if UNION_H5AD.exists() and not args.force:
        print(f"Cached union h5ad present: {UNION_H5AD} — skipping "
              f"(pass --force to rebuild)")
        return

    import anndata as ad
    import scanpy as sc

    print("Step 03 — building CM finetune-data union")
    print("\n[1/4] per-dataset gene maps")
    maps = build_gene_maps()

    print("\n[2/4] per-dataset sample extraction + drug resolution")
    compounds = []
    for ds in CM_DATASETS:
        meta, counts, gene_ids = _datasets.extract_samples(ds)
        agg = meta.groupby(["compound_name", "compound_code"],
                           dropna=False).size().reset_index(name="n_samples")
        for _, r in agg.iterrows():
            smi = meta[(meta["compound_name"] == r["compound_name"])
                       & (meta["compound_code"] == r["compound_code"])
                       ]["smiles"].dropna()
            smi = smi.iloc[0] if len(smi) else None
            compounds.append((ds, r["compound_name"], r["compound_code"],
                              smi, int(r["n_samples"])))
        print(f"  {ds}: {meta.shape[0]} samples, {len(agg)} compounds")
    manifest = build_drug_manifest(compounds)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(MANIFEST, index=False)
    drug_idx_of = {}
    for _, r in manifest.iterrows():
        drug_idx_of[(r["dataset"], str(r["compound_code"]) or r["compound_name"])] \
            = int(r["drug_idx"])

    print("\n[3/4] assemble union h5ad")
    parts = []
    for ds in CM_DATASETS:
        obs, X, X_ctl = process_dataset(ds, maps[ds], drug_idx_of)
        parts.append((ds, obs, X, X_ctl))
        print(f"  {ds}: {len(obs)} treated samples (HG-resolved)")

    # schema: var = 978 landmarks in L1000 gene_idx order
    meta_genes = gene_mapping.load_l1000_meta()
    var = pd.DataFrame({
        "gene_id": meta_genes["gene_id"].values.astype(str),
        "gene_symbol": meta_genes["gene_symbol"].values,
    }, index=meta_genes["gene_id"].values.astype(str))

    obs_all = pd.concat([p[1] for p in parts], ignore_index=True)
    X_all = np.concatenate([p[2] for p in parts], axis=0)
    X_ctl_all = np.concatenate([p[3] for p in parts], axis=0)

    adata = ad.AnnData(X_all.astype(np.float32), obs=obs_all, var=var)
    adata.obsm["X_ctl"] = X_ctl_all.astype(np.float32)

    # provenance / summary
    summary = obs_all.groupby("dataset_id").agg(n_samples=("well", "size"),
                                                n_drugs=("pert_idx", "nunique"))
    summary.to_csv(SUMMARY)
    print("\nDataset summary -> " + str(SUMMARY))
    print(summary)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(UNION_H5AD)
    print(f"\nUnion h5ad: {UNION_H5AD} ({adata.n_obs} samples x {adata.n_vars} genes)")

    print("\n[4/4] schema gate")
    if not args.skip_schema_gate:
        import schema_gate
        schema_gate.run(UNION_H5AD)
    else:
        print("  skipped (--skip-schema-gate)")
    print("Step 03 complete.")


if __name__ == "__main__":
    main()
