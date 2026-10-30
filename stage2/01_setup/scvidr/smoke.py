#!/usr/bin/env python
"""scVIDR setup smoke: 2,000-gene synthetic train + continuous-dose prediction."""
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
MODEL = ROOT / "pert_model" / "scVIDR" / "source"
OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(MODEL / "vidr"))

from utils import prepare_cont_data  # noqa: E402
from vidr import VIDR  # noqa: E402


def main():
    rng = np.random.default_rng(0)
    n_genes = 2000
    doses = [0, 1, 10, 100]
    cell_types = ["A", "B", "C"]
    xs, obs = [], []
    for cell in cell_types:
        for dose in doses:
            n = 20
            xs.append(rng.poisson(1.0, size=(n, n_genes)).astype("float32"))
            obs.append(pd.DataFrame(
                {"celltype": cell, "Dose": float(dose), "dose": str(dose)},
                index=[f"{cell}_{dose}_{i}" for i in range(n)],
            ))
    adata = ad.AnnData(X=np.vstack(xs), obs=pd.concat(obs))
    adata.var_names = [f"G{i}" for i in range(n_genes)]

    train, test = prepare_cont_data(
        adata, "celltype", "dose", "Dose", "C", 0.0, normalized=False
    )
    model = VIDR(train, linear_decoder=False)
    model.train(max_epochs=2, batch_size=64, early_stopping=False)
    pred, delta = model.predict(
        ctrl_key="0",
        treat_key="100",
        cell_type_to_predict="C",
        continuous=True,
        doses=doses,
        low_dose=False,
    )
    np.savez(
        OUT / "pred_doses.npz",
        **{str(dose): arr.X for dose, arr in pred.items()},
    )
    meta = {
        "model": "scVIDR",
        "n_genes": n_genes,
        "train_shape": list(train.shape),
        "test_shape": list(test.shape),
        "predicted_doses": sorted(pred.keys()),
        "predicted_shapes": {str(k): list(v.shape) for k, v in pred.items()},
        "note": "continuous dose is a log-linear latent interpolation, not a learned curve",
    }
    (OUT / "meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
