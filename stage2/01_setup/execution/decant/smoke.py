#!/usr/bin/env python
"""DECANT setup smoke: 2,000-gene mechanism model with two dose values.

The repository ships no processed cell bags or checkpoints, so this is a
synthetic forward pass that proves the model and dose context run.
"""

import sys

# Keep `pert_model/*/source/` byte-identical to its download (ADR 0030): this
# smoke imports model code from the vendored tree, so byte-code must not land
# there. Set before any other import; the run.sh wrapper exports the same flag.
sys.dont_write_bytecode = True

import importlib.util
import json
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent


def find_repo_root(start: Path) -> Path:
    """Walk up to the directory containing .git (ADR 0030 §4)."""
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RuntimeError(f"cannot find the repository root above {start}")


ROOT = find_repo_root(HERE)
MODEL = ROOT / "pert_model" / "DECANT" / "source"
OUT = HERE / "output"
OUT.mkdir(parents=True, exist_ok=True)
SCRIPT = MODEL / "decant" / "scripts" / "train_decant_mechanism_model.py"


def load_module():
    spec = importlib.util.spec_from_file_location("decant_model", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["decant_model"] = module
    spec.loader.exec_module(module)
    return module


def main():
    torch.manual_seed(0)
    m = load_module()
    model = m.Phase1MechanismTopologyModel(
        n_genes=2000,
        n_programs=50,
        n_cell_lines=3,
        token_dim=32,
        gene_hidden=32,
        prog_hidden=32,
        pool_hidden=32,
        drug_dim=16,
        bag_dim=16,
        ctx_dim=8,
        decoder_hidden=32,
        dropout=0.1,
    )
    model.eval()
    b, n = 4, 8
    tg = torch.randn(b, n, 2000)
    tp = torch.randn(b, n, 50)
    tm = torch.ones(b, n, dtype=torch.bool)
    cg = torch.randn(b, n, 2000)
    cp = torch.randn(b, n, 50)
    cm = torch.ones(b, n, dtype=torch.bool)
    cell_idx = torch.tensor([0, 1, 2, 0])
    time = torch.zeros(b)
    with torch.no_grad():
        out1 = model.forward_view(
            tg, tp, tm, cg, cp, cm, cell_idx,
            torch.log1p(torch.full((b,), 1.0)), time,
        )
        out2 = model.forward_view(
            tg, tp, tm, cg, cp, cm, cell_idx,
            torch.log1p(torch.full((b,), 100.0)), time,
        )
    p1 = out1["pred_gene"].numpy()
    p2 = out2["pred_gene"].numpy()
    np.save(OUT / "pred_gene_dose1.npy", p1)
    np.save(OUT / "pred_gene_dose100.npy", p2)
    meta = {
        "model": "DECANT",
        "n_genes": 2000,
        "params": sum(p.numel() for p in model.parameters()),
        "pred_shape": list(p1.shape),
        "dose_effect_l1": float(np.abs(p1 - p2).mean()),
        "note": "synthetic forward; no processed cell bags or checkpoint published",
    }
    (OUT / "meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
