#!/usr/bin/env python
"""biolord setup smoke: train a small model on GSM_new and predict.

GSM_new is a chemical/dose screen (5 drugs, 7 doses, A549). biolord has no
public checkpoint, so this is a short from-scratch run, not a performance test.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import torch

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(parents=True, exist_ok=True)
DATA = ROOT / "pert_model/biolord/source/data/GSM_new.h5ad"

sys.path.insert(0, str(ROOT / "pert_model/biolord/source/src"))
import biolord  # noqa: E402

print("[biolord] loading", DATA)
adata = sc.read_h5ad(DATA)
adata.obs["drug"] = adata.obs["drug"].astype("category")
adata.obs["dose_val"] = adata.obs["dose_val"].astype(float)
print("[biolord] shape", adata.shape)
print("[biolord] drugs", adata.obs["drug"].value_counts().to_dict())
print("[biolord] doses", sorted(adata.obs["dose_val"].unique().tolist()))

biolord.Biolord.setup_anndata(
    adata,
    ordered_attributes_keys=["dose_val"],
    categorical_attributes_keys=["drug"],
)

module_params = {
    "decoder_width": 256,
    "decoder_depth": 2,
    "attribute_nn_width": 128,
    "attribute_nn_depth": 2,
    "n_latent_attribute_categorical": 4,
    "gene_likelihood": "normal",
    "reconstruction_penalty": 1e2,
    "unknown_attribute_penalty": 1e1,
    "unknown_attribute_noise_param": 1e-1,
    "attribute_dropout_rate": 0.1,
    "use_batch_norm": False,
    "use_layer_norm": False,
    "seed": 42,
}

model = biolord.Biolord(
    adata=adata,
    n_latent=16,
    model_name="gsm_smoke",
    module_params=module_params,
    train_classifiers=False,
    split_key="split",
)

use_gpu = torch.cuda.is_available()
print("[biolord] training on", "gpu" if use_gpu else "cpu")
model.train(
    max_epochs=3,
    batch_size=512,
    accelerator="gpu" if use_gpu else "cpu",
    device=1 if use_gpu else "cpu",
    plan_kwargs={
        "n_epochs_warmup": 0,
        "latent_lr": 1e-3,
        "latent_wd": 1e-4,
        "decoder_lr": 1e-3,
        "decoder_wd": 1e-4,
        "attribute_nn_lr": 1e-2,
        "attribute_nn_wd": 4e-8,
        "step_size_lr": 45,
        "cosine_scheduler": False,
    },
    early_stopping=False,
    enable_checkpointing=False,
    num_workers=0,
)

# Counterfactual predictions from control cells in the training split.
ctrl = adata[(adata.obs["split"] == "train") & (adata.obs["control"] == 1)].copy()
ctrl = ctrl[:128].copy()
print("[biolord] source control cells", ctrl.shape)
# biolord's counterfactual API only accepts categorical target attributes;
# dose_val is ordered and is therefore held at each source cell's value.
pred = model.compute_prediction_adata(
    adata,
    ctrl,
    target_attributes=["drug"],
)
X = pred.X.toarray() if hasattr(pred.X, "toarray") else np.asarray(pred.X)
print("[biolord] prediction", X.shape, "range", float(X.min()), float(X.max()))
np.save(OUT / "pred_mean.npy", X)
pred.obs.to_csv(OUT / "pred_obs.csv")

model.save(str(OUT / "model"), overwrite=True, save_anndata=False)
print("[biolord] wrote", OUT / "pred_mean.npy")
