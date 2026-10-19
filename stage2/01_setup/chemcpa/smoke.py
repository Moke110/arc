#!/usr/bin/env python
"""chemCPA setup smoke: load a public checkpoint and run one forward pass.

Finding (2026-10-18): every checkpoint in the public `chemCPA_models.zip` is a
raw `ComPert` state dict trained on the **L1000 977-gene** space (or its 2,000-
gene extension), with 188 drugs and 3 covariates. None is a 5,000-gene or
20,570-gene model. The downloaded `combo_sciplex` data is 5,000 genes, so a real
chemical prediction on it is not possible without retraining.

This script therefore does the documented C2 fallback: build the exact ComPert
architecture, load the state dict, and run a forward pass on a synthetic
L1000-shaped batch. It proves checkpoint/architecture compatibility; it is not a
biological prediction.
"""
from pathlib import Path
import json
import sys

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
OUT.mkdir(parents=True, exist_ok=True)
WEIGHTS = ROOT / "pert_model/chemcpa/weights/chemCPA_models"
DEFAULT_CKPT = WEIGHTS / "27b401db1845eea26c102fb614df9c33.pt"

sys.path.insert(0, str(ROOT / "pert_model/chemcpa"))
from chemCPA.model import ComPert  # noqa: E402


def main():
    ckpt_path = Path(__import__("os").environ.get("CHEMCPA_CKPT", DEFAULT_CKPT))
    print("[chemcpa] checkpoint", ckpt_path)
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    if not isinstance(ckpt, tuple):
        raise RuntimeError(f"unexpected checkpoint type: {type(ckpt)}")
    state_dict, cfg = ckpt[0], ckpt[3]
    emb = state_dict["drug_embeddings.weight"]
    print("[chemcpa] num_genes", cfg["num_genes"], "num_drugs", cfg["num_drugs"],
          "num_covariates", cfg["num_covariates"], "drug_emb", tuple(emb.shape))

    device = "cuda" if torch.cuda.is_available() else "cpu"
    drug_embeddings = torch.nn.Embedding(cfg["num_drugs"], emb.shape[1])
    # The public state dicts carry no covariate-embedding tensors even though
    # their saved config says [3]; they were trained without cell covariates.
    num_covariates = cfg["num_covariates"] if any(
        k.startswith("covariates_embeddings") for k in state_dict
    ) else [0]
    model = ComPert(
        num_genes=cfg["num_genes"],
        num_drugs=cfg["num_drugs"],
        num_covariates=num_covariates,
        device=device,
        seed=cfg["seed"],
        patience=cfg["patience"],
        doser_type=cfg["doser_type"],
        decoder_activation=cfg["decoder_activation"],
        hparams=dict(cfg["hparams"]),
        drug_embeddings=drug_embeddings,
        use_drugs_idx=cfg["use_drugs_idx"],
        append_layer_width=None,
    )
    missing, unexpected = model.load_state_dict(state_dict, strict=False)
    print("[chemcpa] load_state_dict missing", len(missing), "unexpected", len(unexpected))
    if missing:
        raise RuntimeError(f"checkpoint is missing model tensors: {missing[:8]}")
    if unexpected:
        raise RuntimeError(f"checkpoint has unexpected tensors: {unexpected[:8]}")
    model.to(device)
    model.eval()

    b = 8
    genes = torch.randn(b, cfg["num_genes"], device=device)
    drugs_idx = torch.randint(0, cfg["num_drugs"], (b, 1), device=device)
    dosages = torch.rand(b, 1, device=device) * 10
    if num_covariates[0] > 0:
        covariates = [
            F.one_hot(
                torch.randint(0, num_covariates[0], (b,), device=device),
                num_classes=num_covariates[0],
            ).float()
        ]
    else:
        covariates = []
    with torch.no_grad():
        recon, cell_drug = model.predict(
            genes=genes, drugs_idx=drugs_idx, dosages=dosages, covariates=covariates
        )
    pred = recon[:, : cfg["num_genes"]].detach().cpu().numpy()
    print("[chemcpa] prediction", pred.shape, "range", float(pred.min()), float(pred.max()))
    np.save(OUT / "pred.npy", pred)
    (OUT / "meta.json").write_text(
        json.dumps(
            {
                "checkpoint": ckpt_path.name,
                "mode": "load-only-forward",
                "num_genes": cfg["num_genes"],
                "num_drugs": cfg["num_drugs"],
                "num_covariates_model": num_covariates,
                "note": "public chemCPA checkpoints are L1000 977/2000-gene models, "
                        "not 5,000/20,570; real-data inference on combo_sciplex would "
                        "require retraining",
            },
            indent=2,
        )
    )
    print("[chemcpa] wrote", OUT / "pred.npy")


if __name__ == "__main__":
    main()
