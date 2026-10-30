#!/usr/bin/env python
"""Inspect the downloaded chemCPA checkpoints without building a model.

Each public checkpoint is a 5-tuple
    (state_dict, optimizer_states, scheduler_states, model_config, history)
where `model_config` carries num_genes / num_drugs / num_covariates / hparams.
This script records the shapes so a checkpoint can be matched to a dataset.
"""
from pathlib import Path
import json
import sys

import torch

ROOT = Path(__file__).resolve().parents[3]
WEIGHTS = ROOT / "pert_model/chemCPA/source/weights/chemCPA_models"
OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "pert_model/chemCPA/source"))

rows = []
for p in sorted(WEIGHTS.glob("*.pt")):
    row = {"file": p.name}
    try:
        obj = torch.load(p, map_location="cpu", weights_only=False)
        if isinstance(obj, tuple):
            state_dict, cfg = obj[0], obj[3]
        else:
            state_dict, cfg = obj.get("state_dict", {}), obj.get("hyper_parameters", {}).get("dataset_config", {})
        emb = state_dict.get("drug_embeddings.weight")
        row.update(
            {
                "status": "ok",
                "num_genes": cfg.get("num_genes"),
                "num_drugs": cfg.get("num_drugs"),
                "num_covariates": cfg.get("num_covariates"),
                "use_drugs_idx": cfg.get("use_drugs_idx"),
                "doser_type": cfg.get("doser_type"),
                "drug_emb_shape": list(emb.shape) if emb is not None else None,
                "has_covariate_embeddings": any(
                    k.startswith("covariates_embeddings") for k in state_dict
                ),
            }
        )
    except Exception as e:  # pragma: no cover - diagnostic
        row.update({"status": "error", "error": f"{type(e).__name__}: {e}"})
    rows.append(row)
    print(json.dumps(row, ensure_ascii=False))

(OUT / "checkpoint_inventory.json").write_text(json.dumps(rows, indent=2))
print("wrote", OUT / "checkpoint_inventory.json")
