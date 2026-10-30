#!/usr/bin/env python
"""CPA setup smoke: load the public Kang checkpoint and predict.

This is deliberately a native-asset smoke, not a benchmark. It uses Kang
(IFN-beta stimulation) because that is the only public CPA checkpoint that
matches a downloaded dataset; it says nothing about cold-drug performance.
"""

import sys

# Keep `pert_model/*/source/` byte-identical to its download (ADR 0030): this
# smoke imports model code from the vendored tree, so byte-code must not land
# there. Set before any other import; the run.sh wrapper exports the same flag.
sys.dont_write_bytecode = True

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent


def find_repo_root(start: Path) -> Path:
    """Walk up to the directory containing .git (ADR 0030 §4)."""
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RuntimeError(f"cannot find the repository root above {start}")


ROOT = find_repo_root(HERE)
OUT = HERE / "output"
OUT.mkdir(parents=True, exist_ok=True)
CKPT = ROOT / "pert_model/CPA/source/weights/Pretrained_Models/Kang"
DATA = ROOT / "pert_model/CPA/source/data/kang_normalized_hvg.h5ad"

# cpa/__init__.py resolves its own version through importlib.metadata, which
# fails for a plain source checkout. The old tree shipped a fake `pydist/`
# dist-info package for this; that stub was deleted in ticket 06, so the
# fallback lives here at runtime and never lands in git.
import importlib.metadata as _metadata  # noqa: E402

_real_version = _metadata.version


def _version(name):
    try:
        return _real_version(name)
    except _metadata.PackageNotFoundError:
        if name == "cpa-tools":
            return "0.0.0"
        raise


_metadata.version = _version
sys.path.insert(0, str(ROOT / "pert_model/CPA/source"))
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
