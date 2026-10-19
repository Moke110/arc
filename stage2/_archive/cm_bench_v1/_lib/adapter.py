"""
The single benchmark seam: the `predict` adapter contract.

Every entrant — model or baseline — implements:

    predict(train_cases, test_cases) -> Δ matrix (n_test, n_genes)

`train_cases` / `test_cases` are `Cases` bundles carrying the per-plate baseline
transcriptome, drug identity (SMILES where available), dose (log10 µM) and — for
train only — the observed Δ target. The adapter owns finetuning *and*
prediction; the workflow never reaches into model internals. `NoChangeAdapter`
returns zeros and doubles as the fake adapter used to exercise the whole
evaluate path on CPU without any model environment.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class Cases:
    baseline: np.ndarray                       # (n, G) per-plate baseline
    drug: list                                 # (n,) drug name
    dataset: list                              # (n,) dataset id
    dose_log: np.ndarray                       # (n,) log10 uM
    group_key: list                            # (n,) (dataset, drug, dose_level)
    smiles: Optional[list] = None              # (n,) SMILES or None
    delta: Optional[np.ndarray] = None         # (n, G) train target only

    def __len__(self):
        return self.baseline.shape[0]

    @property
    def n_genes(self):
        return self.baseline.shape[1]

    @classmethod
    def from_frame(cls, meta, baseline, delta=None):
        keys = list(zip(meta["dataset"], meta["drug"], meta["dose_level"]))
        return cls(
            baseline=np.asarray(baseline, dtype=float),
            drug=list(meta["drug"]),
            dataset=list(meta["dataset"]),
            dose_log=meta["log10_dose"].to_numpy(dtype=float),
            group_key=keys,
            smiles=list(meta["smiles"]) if "smiles" in meta else None,
            delta=None if delta is None else np.asarray(delta, dtype=float),
        )


class NoChangeAdapter:
    """Δ = 0. Null baseline and the CPU fake adapter for evaluation tests."""

    name = "no-change"
    supports_dose = True

    def predict(self, train_cases, test_cases):
        return np.zeros((len(test_cases), test_cases.n_genes), dtype=float)
