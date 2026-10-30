#!/usr/bin/env python
"""Inspect the downloaded chemCPA checkpoints without building a model.

Each public checkpoint is a 5-tuple
    (state_dict, optimizer_states, scheduler_states, model_config, history)
where `model_config` carries num_genes / num_drugs / num_covariates / hparams.
This script records the shapes so a checkpoint can be matched to a dataset.
"""

import sys

# Keep `pert_model/*/source/` byte-identical to its download (ADR 0030): this
# smoke imports model code from the vendored tree, so byte-code must not land
# there. Set before any other import; the run.sh wrapper exports the same flag.
sys.dont_write_bytecode = True

from pathlib import Path
import json

import torch

HERE = Path(__file__).resolve().parent


def find_repo_root(start: Path) -> Path:
    """Walk up to the directory containing .git (ADR 0030 §4)."""
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RuntimeError(f"cannot find the repository root above {start}")


ROOT = find_repo_root(HERE)
WEIGHTS = ROOT / "pert_model/chemCPA/source/weights/chemCPA_models"
OUT = HERE / "output"
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
