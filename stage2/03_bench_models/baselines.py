"""Mandatory non-neural baselines for the Stage 2 benchmark (ticket 12).

All baselines are fit on the training wells of one fold and predict a per-well
Δ.  They are invoked by `adapters/baselines/run.py` through the same
subprocess seam as the real models; the runner averages their per-well outputs
within `group_id` before scoring.

Baselines
---------
- ``delta_mean``: the training group-mean Δ, identical for every test group.
- ``delta_mean_dose``: the training group-mean Δ per dose.
- ``ridge_fp``: ridge regression from an ECFP4 fingerprint + log-dose to the
  training group-mean Δ.
- ``nearest_drug``: the group-mean Δ of the nearest training drug (Tanimoto on
  ECFP4), at the same dose when available.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import Ridge

from cases import save_npz


def morgan_fp(smiles: str, n_bits: int = 2048, radius: int = 2) -> np.ndarray:
    """ECFP4-style Morgan fingerprint as a float32 bit vector."""

    arr = np.zeros(n_bits, dtype=np.float32)
    mol = Chem.MolFromSmiles(str(smiles)) if str(smiles) else None
    if mol is None:
        return arr
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = generator.GetFingerprint(mol)
    for bit in fp.GetOnBits():
        arr[int(bit)] = 1.0
    return arr


def _group_table(cases: np.lib.npyio.NpzFile) -> dict:
    """Collapse per-well train/valid rows into training group means."""

    group_ids = np.asarray(cases["group_id"]).astype(str)
    unique = np.unique(group_ids)
    delta = np.asarray(cases["delta"], dtype=np.float32)
    if delta.shape[1] != 18164:
        raise RuntimeError(f"unexpected Δ width {delta.shape[1]}")
    group_delta = np.empty((len(unique), delta.shape[1]), dtype=np.float64)
    smiles = []
    doses = []
    drugs = []
    for i, group in enumerate(unique):
        rows = np.flatnonzero(group_ids == group)
        group_delta[i] = delta[rows].mean(axis=0, dtype=np.float64)
        smiles.append(str(cases["smiles"][rows[0]]))
        doses.append(float(cases["dose_um"][rows[0]]))
        drugs.append(str(cases["drug_key"][rows[0]]))
    return {
        "group_ids": unique,
        "group_delta": group_delta,
        "smiles": np.asarray(smiles, dtype=str),
        "doses": np.asarray(doses, dtype=np.float64),
        "drugs": np.asarray(drugs, dtype=str),
    }


def _ridge_features(
    smiles: np.ndarray, doses: np.ndarray, n_bits: int, radius: int
) -> np.ndarray:
    fps = np.stack([morgan_fp(s, n_bits, radius) for s in smiles])
    dose_feature = np.log10(doses.astype(np.float64) + 1.0)[:, None]
    return np.concatenate([fps.astype(np.float64), dose_feature], axis=1)


def _nearest_group(
    smiles: str, dose: float, train: dict, n_bits: int, radius: int
) -> np.ndarray:
    query = morgan_fp(smiles, n_bits, radius)
    train_fp = train["fp"]
    intersection = train_fp @ query
    union = train_fp.sum(axis=1) + query.sum() - intersection
    similarity = np.divide(
        intersection, union, out=np.zeros_like(intersection), where=union > 0
    )
    nearest = int(np.argmax(similarity))
    nearest_drug = str(train["drugs"][nearest])
    rows = np.flatnonzero(train["group_drug"] == nearest_drug)
    if not len(rows):
        return train["fallback"]
    doses = train["group_dose"][rows]
    pick = rows[int(np.argmin(np.abs(doses - float(dose))))]
    return train["group_delta"][pick]


def fit_baseline(
    name: str,
    train_path: Path,
    valid_path: Path,
    config: dict,
    run_dir: Path,
) -> Path:
    """Fit one baseline and write a checkpoint directory."""

    run_dir = Path(run_dir)
    checkpoint = run_dir / "checkpoint"
    checkpoint.mkdir(parents=True, exist_ok=True)
    params = config.get("baseline_params", {})
    n_bits = int(params.get("fp_bits", 2048))
    radius = int(params.get("fp_radius", 2))
    alpha = float(params.get("ridge_alpha", 1.0))

    with np.load(train_path, allow_pickle=False) as cases:
        table = _group_table(cases)
    group_delta = table["group_delta"]
    meta = {
        "baseline": name,
        "n_groups": int(len(table["group_ids"])),
        "n_genes": int(group_delta.shape[1]),
        "fp_bits": n_bits,
        "fp_radius": radius,
        "ridge_alpha": alpha,
    }

    if name == "delta_mean":
        save_npz(checkpoint / "params.npz", {"mean": group_delta.mean(axis=0).astype(np.float32)})
    elif name == "delta_mean_dose":
        doses = np.unique(table["doses"])
        means = np.stack(
            [group_delta[table["doses"] == d].mean(axis=0) for d in doses]
        ).astype(np.float32)
        save_npz(checkpoint / "params.npz", {"doses": doses.astype(np.float64), "means": means})
        meta["doses"] = doses.tolist()
    elif name == "ridge_fp":
        X = _ridge_features(table["smiles"], table["doses"], n_bits, radius)
        model = Ridge(alpha=alpha, fit_intercept=True, solver="cholesky")
        model.fit(X, group_delta)
        save_npz(
            checkpoint / "params.npz",
            {
                "coef": np.asarray(model.coef_, dtype=np.float32),
                "intercept": np.asarray(model.intercept_, dtype=np.float32),
            },
        )
        meta["n_features"] = int(X.shape[1])
    elif name == "nearest_drug":
        unique_drugs = np.unique(table["drugs"])
        drug_fp = np.stack(
            [
                morgan_fp(table["smiles"][np.flatnonzero(table["drugs"] == d)[0]], n_bits, radius)
                for d in unique_drugs
            ]
        )
        save_npz(
            checkpoint / "params.npz",
            {
                "drugs": unique_drugs.astype(str),
                "fp": drug_fp.astype(np.float32),
                "group_drug": table["drugs"].astype(str),
                "group_dose": table["doses"].astype(np.float64),
                "group_delta": group_delta.astype(np.float32),
                "fallback": group_delta.mean(axis=0).astype(np.float32),
            },
        )
        meta["n_train_drugs"] = int(len(unique_drugs))
    else:
        raise ValueError(f"unknown baseline {name}")

    (checkpoint / "meta.json").write_text(
        json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (run_dir / "checkpoint_path.txt").write_text(str(checkpoint.resolve()) + "\n", encoding="utf-8")
    return checkpoint


def predict_baseline(
    name: str,
    cases_path: Path,
    checkpoint: Path,
    out_path: Path,
) -> np.ndarray:
    """Predict one Δ row per case (well), in the bundle's row order."""

    checkpoint = Path(checkpoint)
    with np.load(checkpoint / "params.npz", allow_pickle=False) as params:
        with np.load(cases_path, allow_pickle=False) as cases:
            group_ids = np.asarray(cases["group_id"]).astype(str)
            unique, inverse = np.unique(group_ids, return_inverse=True)
            n_genes = int(params["mean"].shape[0]) if "mean" in params.files else 18164
            pred_groups = np.zeros((len(unique), n_genes), dtype=np.float32)

            if name == "delta_mean":
                pred_groups[:] = params["mean"]
            elif name == "delta_mean_dose":
                doses = np.asarray(params["doses"], dtype=np.float64)
                means = np.asarray(params["means"], dtype=np.float32)
                for i, group in enumerate(unique):
                    rows = np.flatnonzero(group_ids == group)
                    dose = float(cases["dose_um"][rows[0]])
                    j = int(np.argmin(np.abs(doses - dose)))
                    pred_groups[i] = means[j]
            elif name == "ridge_fp":
                coef = np.asarray(params["coef"], dtype=np.float32)
                intercept = np.asarray(params["intercept"], dtype=np.float32)
                n_bits = int(coef.shape[1] - 1)
                for i, group in enumerate(unique):
                    rows = np.flatnonzero(group_ids == group)
                    smiles = str(cases["smiles"][rows[0]])
                    dose = float(cases["dose_um"][rows[0]])
                    fp = morgan_fp(smiles, n_bits, 2)
                    x = np.concatenate([fp, [np.log10(dose + 1.0)]]).astype(np.float32)
                    pred_groups[i] = coef @ x + intercept
            elif name == "nearest_drug":
                train = {
                    "drugs": np.asarray(params["drugs"]).astype(str),
                    "fp": np.asarray(params["fp"], dtype=np.float32),
                    "group_drug": np.asarray(params["group_drug"]).astype(str),
                    "group_dose": np.asarray(params["group_dose"], dtype=np.float64),
                    "group_delta": np.asarray(params["group_delta"], dtype=np.float32),
                    "fallback": np.asarray(params["fallback"], dtype=np.float32),
                }
                for i, group in enumerate(unique):
                    rows = np.flatnonzero(group_ids == group)
                    smiles = str(cases["smiles"][rows[0]])
                    dose = float(cases["dose_um"][rows[0]])
                    pred_groups[i] = _nearest_group(smiles, dose, train, 2048, 2)
            else:
                raise ValueError(f"unknown baseline {name}")

            pred = pred_groups[inverse]
    np.save(out_path, np.ascontiguousarray(pred, dtype=np.float32))
    return pred
