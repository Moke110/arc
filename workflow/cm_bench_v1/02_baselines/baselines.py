"""
Mandatory baselines (tickets 05) implementing the `predict` adapter seam.

  no-change     Δ = 0                        (null; also the CPU fake adapter)
  mean-effect   Δ = mean train Δ             (global)
  ridge-fp      Morgan/ECFP fingerprint + log-dose -> ridge -> Δ
  nearest-drug  copy the mean Δ of the most fingerprint-similar train drug

The fingerprint is Morgan radius 2 with 2048 bits; dose is log10 µM. Baselines
that need chemistry fall back to the global mean when a SMILES is missing.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from adapter import NoChangeAdapter  # noqa: E402,F401  (re-exported)

FP_RADIUS = 2
FP_BITS = 2048


def _fp(smiles):
    if not smiles:
        return None
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, FP_RADIUS, nBits=FP_BITS)
    arr = np.zeros(FP_BITS, dtype=np.float64)
    for bit in fp.GetOnBits():
        arr[bit] = 1.0
    return arr


def _tanimoto(a, b):
    inter = float((a * b).sum())
    union = float(a.sum() + b.sum() - inter)
    return inter / union if union > 0 else 0.0


def _train_mean(train):
    if train.delta is None or len(train) == 0:
        return np.zeros(train.n_genes)
    return train.delta.mean(axis=0)


class MeanEffectBaseline:
    name = "mean-effect"
    supports_dose = False

    def predict(self, train, test):
        mean = _train_mean(train)
        return np.tile(mean, (len(test), 1))


class RidgeFingerprintBaseline:
    name = "ridge-fp"
    supports_dose = True

    def __init__(self, alpha=1.0):
        self.alpha = alpha

    def _features(self, cases):
        X = np.zeros((len(cases), FP_BITS + 1))
        for i, (smi, dose) in enumerate(zip(cases.smiles or [None] * len(cases),
                                            cases.dose_log)):
            fp = _fp(smi)
            if fp is not None:
                X[i, :FP_BITS] = fp
            X[i, FP_BITS] = dose
        return X

    def predict(self, train, test):
        from sklearn.linear_model import Ridge

        mean = _train_mean(train)
        Xtr = self._features(train)
        if train.delta is None or len(train) == 0 or not Xtr.any():
            return np.tile(mean, (len(test), 1))
        model = Ridge(alpha=self.alpha)
        model.fit(Xtr, train.delta)
        pred = model.predict(self._features(test))
        return pred


class NearestDrugBaseline:
    name = "nearest-drug"
    supports_dose = False

    def predict(self, train, test):
        mean = _train_mean(train)
        # mean Δ per train drug (fingerprint + delta)
        by_drug = {}
        for i in range(len(train)):
            d = train.drug[i]
            if d not in by_drug:
                by_drug[d] = {"fp": _fp((train.smiles or [None] * len(train))[i]),
                              "sum": np.zeros(train.n_genes), "n": 0}
            by_drug[d]["sum"] += train.delta[i]
            by_drug[d]["n"] += 1
        drug_fp = {d: v["fp"] for d, v in by_drug.items()}
        drug_mean = {d: v["sum"] / v["n"] for d, v in by_drug.items()}
        pred = np.tile(mean, (len(test), 1))
        for i, smi in enumerate(test.smiles or [None] * len(test)):
            q = _fp(smi)
            if q is None:
                continue
            best, best_sim = None, -1.0
            for d, fp in drug_fp.items():
                if fp is None:
                    continue
                sim = _tanimoto(q, fp)
                if sim > best_sim:
                    best, best_sim = d, sim
            if best is not None:
                pred[i] = drug_mean[best]
        return pred


BASELINES = {
    "no-change": NoChangeAdapter,
    "mean-effect": MeanEffectBaseline,
    "ridge-fp": RidgeFingerprintBaseline,
    "nearest-drug": NearestDrugBaseline,
}
