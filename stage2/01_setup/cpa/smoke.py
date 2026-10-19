#!/usr/bin/env python
"""CPA setup smoke: load the public Kang checkpoint and predict.

This is deliberately a native-asset smoke, not a benchmark. It uses Kang
(IFN-beta stimulation) because that is the only public CPA checkpoint that
matches a downloaded dataset; it says nothing about cold-drug performance.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(parents=True, exist_ok=True)
CKPT = ROOT / "pert_model/cpa/weights/Pretrained_Models/Kang"
DATA = ROOT / "pert_model/cpa/data/kang_normalized_hvg.h5ad"

# cpa/__init__.py reads its own package version; we run from a source checkout
# rather than an installed wheel, so provide a tiny local dist-info.
sys.path.insert(0, str(Path(__file__).resolve().parent / "pydist"))
sys.path.insert(0, str(ROOT / "pert_model/cpa"))
import cpa  # noqa: E402

print("[cpa] loading", DATA)
adata = sc.read_h5ad(DATA)
print("[cpa] raw", adata.shape)

# A small, cell-type-balanced subset keeps the CPU smoke quick while retaining
# every category the checkpoint knows about.
keep = (
    adata.obs.groupby("cell_type", observed=True)
    .head(250)
    .index
)
adata = adata[keep].copy()
adata.obs["dose"] = 1.0
adata.obs["condition"] = adata.obs["condition"].astype("category")
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
print("[cpa] smoke", adata.shape, adata.obs["cell_type"].value_counts().to_dict())

cpa.CPA.setup_anndata(
    adata,
    perturbation_key="condition",
    control_group="ctrl",
    dosage_key="dose",
    categorical_covariate_keys=["cell_type"],
    is_count_data=False,
    max_comb_len=1,
)

print("[cpa] loading checkpoint", CKPT)
model = cpa.CPA.load(dir_path=str(CKPT), adata=adata, use_gpu=False)
model.predict(adata, batch_size=256, n_samples=5, return_mean=True)

if "CPA_pred" not in adata.obsm:
    raise RuntimeError(f"CPA.predict did not write adata.obsm['CPA_pred']; keys={list(adata.obsm)}")

pred = np.asarray(adata.obsm["CPA_pred"])
print("[cpa] prediction", pred.shape, "range", float(pred.min()), float(pred.max()))
np.save(OUT / "cpa_pred.npy", pred)
adata.obs.to_csv(OUT / "cpa_pred_obs.csv")
print("[cpa] wrote", OUT / "cpa_pred.npy")
