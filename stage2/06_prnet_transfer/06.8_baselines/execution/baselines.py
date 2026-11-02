#!/usr/bin/env python
"""Ticket #8 — five-fold simple baselines on the frozen 06.3 input.

The baselines exist so the formal models (tickets #5–#7, scored by #9) can be
paired against them on the **identical main test population**, with the
**identical evaluator** and the **identical truth Δ**.  This module therefore
does not re-implement the metric: it writes prediction directories in the
``06.4_pipeline/execution/predict.py`` format (``pred_expr.npy`` / ``delta.npy``
/ ``meta.csv`` / ``manifest.json``) and scores each one by calling
``ptlib.evaluate.run_evaluation`` — the single evaluator of the study.

Definitions (frozen before any test score is seen; see ``README.md``):

``no_change``
    Δ = 0 for every well, i.e. predicted expression = the fold's fixed consensus
    ``c*`` — the constant-prediction boundary the spec asks for.
``delta_mean``
    The historical training-mean constant baseline.  Pooling: within-condition
    replicate aggregation **first** (condition = ``drug_key|dose_µM``), then an
    equal-weight arithmetic mean over the fold's **training** conditions.  A test
    condition whose exact key also occurs in training uses that training
    condition's mean Δ; otherwise the pooled training mean is the documented
    fallback.  Under the frozen drug-key split a test drug is never a training
    drug, so the fallback is always taken; the run counts and records the
    collisions instead of assuming zero.
``delta_mean_dose``
    Auxiliary: the training condition means pooled per dose (:g of µM); a dose
    absent from training falls back to the pooled training mean.
``ridge_fp``
    Auxiliary: L2 ridge (``alpha=1.0``, fixed a priori, never tuned) from the
    study's own dose encoding — 1024-bit rFCFP4 fingerprint × log10(dose_µM+1),
    exactly the model's ``drug_dim`` label — to the training condition-mean Δ.
    Fitted per fold with ``fit_intercept=True``.

Leakage rule: only the fold's **training** treated responses enter a fit.  The
matched-control means used to form the training Δ are controls, not held-out
responses: they are never split by drug and the scored test truth uses the same
per-plate DMSO definition.  Validation/test treated rows are never read here.

Memory (ticket #8; ``docs/agents/execution-environment.md``): the dense CM log1p
layer is read and aggregated one ``GENE_BLOCK`` of genes at a time, the matched
controls are kept as ``n_plates × n_genes`` plate means, and the two prediction
``.npy`` files are written through ``np.lib.format`` memmaps in ``ROW_BLOCK`` row
blocks.  No ``n_observations × 18_164`` array is ever materialised, and there is
no worker pool (single process).  ``tests/test_memory.py`` pins the block
equivalence, and ``make mem-check`` measures the peak on a row-replicated fixture.

Usage::

    PY=/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
    $PY execution/baselines.py --check
    $PY execution/baselines.py --all \\
        --data-dir /home/chang/arc/stage2/06_prnet_transfer/06.3_input/output \\
        --out-root /home/chang/arc/stage2/06_prnet_transfer

Exit codes: 0 success, 1 contract/data error.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

__all__ = [
    "BASELINES",
    "PRIMARY_BASELINES",
    "AUXILIARY_BASELINES",
    "RIDGE_ALPHA",
    "GENE_BLOCK",
    "ROW_BLOCK",
    "FitTable",
    "PredictionResult",
    "find_repo_root",
    "iter_layer_blocks",
    "plate_control_matrix",
    "build_fit_table",
    "condition_context",
    "ridge_features",
    "fit_baseline",
    "condition_predictions",
    "predict_baseline",
    "write_prediction_dir",
    "score_prediction_dir",
    "run_fold",
    "run_all",
    "check_input",
    "summarise_folds",
    "headline_metrics",
    "arrays_sha256",
    "SCHEMA",
]

SCHEMA = "stage2/06.8_baselines/v1"
BASELINES = ("no_change", "delta_mean", "delta_mean_dose", "ridge_fp")
#: Baselines the spec requires in the main comparison.
PRIMARY_BASELINES = ("no_change", "delta_mean")
#: Baselines the spec suggests as auxiliary; never part of the main comparison.
AUXILIARY_BASELINES = ("delta_mean_dose", "ridge_fp")
#: Fixed a priori (the 03.1 ``ridge_alpha`` default); never tuned on any split.
RIDGE_ALPHA = 1.0
#: The step container whose ``06.8_baselines/output/`` this ticket owns.
STEP_REL = Path("stage2/06_prnet_transfer")

#: Genes per block when the dense CM log1p layer is read/aggregated.  The only
#: large temporary is ``(n_train_rows × GENE_BLOCK)`` float64, never
#: ``n_train_rows × 18_164`` (ticket #8 memory rule; execution-environment.md).
GENE_BLOCK = 1024
#: Scored rows per block when the prediction arrays are written.  Both ``.npy``
#: files are created as ``np.lib.format`` memmaps and filled block by block, so
#: no ``n_test_rows × 18_164`` array is ever materialised.
ROW_BLOCK = 64


# --------------------------------------------------------------------------- #
# repository / pipeline locations (walk to .git; never count parents)
# --------------------------------------------------------------------------- #
class BaselineError(RuntimeError):
    """A contract or data invariant of ticket #8 does not hold."""


def find_repo_root(start: Path | str | None = None) -> Path:
    path = Path(start).resolve() if start is not None else Path(__file__).resolve()
    if path.is_file():
        path = path.parent
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise BaselineError(f"no .git directory at or above {path}")


ROOT = find_repo_root(Path(__file__).resolve())
#: The library of ticket #4: the single training/evaluation code path.
PIPELINE_EXECUTION = ROOT / STEP_REL / "06.4_pipeline" / "execution"
if str(PIPELINE_EXECUTION) not in sys.path:  # pragma: no cover - import plumbing
    sys.path.insert(0, str(PIPELINE_EXECUTION))

from ptlib import PIPELINE_VERSION  # noqa: E402
from ptlib import data as data_mod  # noqa: E402
from ptlib import evaluate as evaluate_mod  # noqa: E402
from ptlib import h5io  # noqa: E402
from ptlib import manifest as manifest_mod  # noqa: E402
from ptlib import metrics as metrics_mod  # noqa: E402


# --------------------------------------------------------------------------- #
# hashing / identity helpers
# --------------------------------------------------------------------------- #
def _sha256_array(array: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(str(contiguous.dtype).encode("ascii"))
    digest.update(str(contiguous.shape).encode("ascii"))
    digest.update(contiguous.tobytes())
    return digest.hexdigest()


def arrays_sha256(arrays: Mapping[str, np.ndarray]) -> str:
    """Content hash of a parameter dict, independent of npz container bytes."""
    digest = hashlib.sha256()
    for key in sorted(arrays):
        digest.update(key.encode("utf-8"))
        digest.update(_sha256_array(np.asarray(arrays[key])).encode("ascii"))
    return digest.hexdigest()


def _sha256_rows(rows: np.ndarray) -> str:
    return _sha256_array(np.asarray(rows, dtype=np.int64))


def _environment() -> dict[str, Any]:
    import sklearn

    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "sklearn": sklearn.__version__,
        "platform": platform.platform(),
    }


# --------------------------------------------------------------------------- #
# fit table: the fold's TRAINING responses only
# --------------------------------------------------------------------------- #
@dataclass
class FitTable:
    """Per-fold training Δ collapsed to drug×dose conditions.

    ``train_rows`` are the fold's ``role=train`` treated observations; the
    per-row Δ is ``treated log1p_cp10k − per-plate matched DMSO mean`` — the same
    control definition the evaluator uses for the scored truth Δ.  Conditions are
    ``drug_key|dose_µM`` (the evaluator's CM condition), and ``condition_delta``
    is the within-condition mean, i.e. replicate aggregation happens **before**
    any pooling.

    Memory (ticket #8): ``condition_delta`` is condition-level
    (``n_conditions × n_genes`` float32) and is the only large array retained;
    the ``n_train_rows × n_genes`` row matrix is never materialised
    (:func:`build_fit_table` aggregates gene block by gene block).
    """

    fold: int
    n_genes: int
    train_rows: np.ndarray
    row_conditions: np.ndarray
    unique_conditions: np.ndarray
    condition_delta: np.ndarray
    condition_drug: np.ndarray
    condition_dose: np.ndarray
    pooled_mean: np.ndarray
    counts: dict[str, Any] = field(default_factory=dict)
    train_rows_sha256: str = ""
    train_delta_sha256: str = ""

    @property
    def n_rows(self) -> int:
        return int(self.train_rows.size)

    @property
    def n_conditions(self) -> int:
        return int(self.unique_conditions.size)

    def condition_lookup(self) -> dict[str, int]:
        return {str(key): i for i, key in enumerate(self.unique_conditions)}


def _row_runs(rows: np.ndarray) -> list[tuple[int, int]]:
    """Split sorted unique *rows* into contiguous ``(start, stop)`` HDF5 ranges.

    Mirrors the row-read contract of ``ptlib.h5io``: a dense frozen matrix is
    read in contiguous ranges, never with one unbounded fancy index.
    """
    runs: list[tuple[int, int]] = []
    start = 0
    for position in range(1, len(rows) + 1):
        if position == len(rows) or rows[position] != rows[position - 1] + 1:
            runs.append((int(rows[start]), int(rows[position - 1]) + 1))
            start = position
    return runs


def _log1p_layer_key(handle: h5py.File, path: Path) -> str:
    for candidate in h5io.LOG1P_KEYS:
        if candidate in handle:
            return candidate
    raise BaselineError(
        f"{path} has no recognised log1p-CP10k layer (looked for {list(h5io.LOG1P_KEYS)})"
    )


def iter_layer_blocks(
    path: Path | str,
    rows: np.ndarray,
    n_genes: int,
    *,
    block_genes: int = GENE_BLOCK,
):
    """Yield ``(g0, g1, block)`` float32 gene blocks of the dense CM log1p layer.

    ``block`` is ``(len(rows), g1 − g0)`` in the given row order.  The full
    ``rows × n_genes`` matrix is never materialised: each block is read straight
    from the contiguous dense dataset and aggregated by the caller.  A CSR layer
    is refused instead of densified — the frozen CM input contract (ticket #3) is
    a dense float32 layer, so an unexpected encoding fails closed rather than
    silently allocating the whole matrix.
    """
    rows = np.asarray(rows, dtype=np.int64)
    if rows.ndim != 1 or rows.size == 0:
        raise BaselineError("layer block read needs a non-empty 1-D row set")
    if np.any(np.diff(rows) <= 0):
        raise BaselineError("layer block read needs sorted, unique row ids")
    path = Path(path)
    with h5py.File(path, "r") as handle:
        key = _log1p_layer_key(handle, path)
        node = handle[key]
        if isinstance(node, h5py.Group):
            raise BaselineError(
                f"{path} stores {key!r} as CSR; the CM input contract is a dense "
                "float32 layer and the baselines read it in gene blocks instead of "
                "densifying the whole matrix"
            )
        shape = tuple(int(x) for x in node.shape)
        if len(shape) != 2 or shape[1] != int(n_genes):
            raise BaselineError(f"{key!r} shape {shape} is not (*, {int(n_genes)})")
        if int(rows.max()) >= shape[0]:
            raise BaselineError(
                f"row {int(rows.max())} is outside {key!r} with {shape[0]} rows"
            )
        runs = _row_runs(rows)
        for g0 in range(0, int(n_genes), int(block_genes)):
            g1 = min(g0 + int(block_genes), int(n_genes))
            pieces = [node[r0:r1, g0:g1] for r0, r1 in runs]
            block = pieces[0] if len(pieces) == 1 else np.vstack(pieces)
            block = np.ascontiguousarray(block, dtype=np.float32)
            if block.shape != (rows.size, g1 - g0):
                raise BaselineError(
                    f"layer block read returned {block.shape}, expected "
                    f"{(rows.size, g1 - g0)}"
                )
            yield g0, g1, block


def plate_control_matrix(
    pool: data_mod.CMControlPool, rows: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Matched-control values compactly: ``(plate_means, plate_index_per_row)``.

    Equivalent to :meth:`ptlib.data.CMControlPool.mean_for_rows` — same
    ``plate_means`` table, same float32 values, same fail-closed missing-plate
    rule — but keeps ``n_plates × n_genes`` instead of materialising
    ``n_rows × n_genes``; the caller indexes it one gene block at a time.
    """
    frozen = pool.frozen
    plate = frozen.cm_obs["plate"].astype(str).to_numpy()
    rows = np.asarray(rows, dtype=np.int64)
    plates = sorted(str(key) for key in pool.plate_means)
    index = {key: i for i, key in enumerate(plates)}
    per_row = np.empty(rows.size, dtype=np.int64)
    for i, row in enumerate(rows):
        key = str(plate[int(row)])
        position = index.get(key)
        if position is None:
            raise BaselineError(
                f"no QC-pass matched DMSO control on plate {key!r} for the fit Δ "
                "(fails closed; one plate's control is never borrowed from another)"
            )
        per_row[i] = position
    matrix = np.stack([pool.plate_means[key] for key in plates]).astype(np.float32)
    if matrix.shape != (len(plates), int(frozen.n_genes)):
        raise BaselineError(
            f"plate control matrix shape {matrix.shape} is not "
            f"({len(plates)}, {frozen.n_genes})"
        )
    return matrix, per_row


def build_fit_table(
    frozen: data_mod.FrozenInput,
    fold_table: data_mod.FoldTable,
    fold: int,
    structures: data_mod.DrugStructures | None = None,
    control_pool: data_mod.CMControlPool | None = None,
    source: data_mod.CMPredictSource | None = None,
) -> FitTable:
    """Read the fold's training treated rows and build the condition fit table.

    Memory-bounded by construction (ticket #8): the dense log1p layer is read and
    aggregated one ``GENE_BLOCK`` at a time, so the peak temporary is
    ``n_train_rows × GENE_BLOCK`` float64 — never ``n_train_rows × n_genes``.
    Fails closed if a row set is empty, if a block read returns an unexpected
    shape/order, or if the condition key disagrees with the (drug, dose) pair.
    """
    structures = structures or frozen.structures()
    source = source or data_mod.CMPredictSource(frozen, fold_table, int(fold), structures)
    train_rows = np.asarray(source.rows_for("train"), dtype=np.int64)
    if train_rows.size == 0:
        raise BaselineError(f"fold {fold}: no training treated rows; a fit is impossible")
    if np.any(np.diff(train_rows) <= 0):
        raise BaselineError(f"fold {fold}: training rows are not sorted and unique")

    keys = np.asarray(source.condition_key(train_rows), dtype=object)
    unique, inverse = np.unique(keys, return_inverse=True)
    inverse = np.asarray(inverse, dtype=np.int64)
    n_conditions = int(unique.size)

    obs = frozen.cm_obs
    drug = obs["drug_key"].astype(str).to_numpy()[train_rows]
    dose = pd.to_numeric(obs["dose_um"], errors="coerce").to_numpy(dtype=np.float64)[train_rows]
    first = np.full(n_conditions, -1, dtype=np.int64)
    for position, condition in enumerate(inverse):
        if first[condition] < 0:
            first[condition] = position
    if np.any(first < 0):
        raise BaselineError(f"fold {fold}: a training condition has no row")
    condition_drug = drug[first]
    condition_dose = dose[first]
    for i, key in enumerate(unique):
        expected = f"{condition_drug[i]}|{condition_dose[i]:g}"
        if str(key) != expected:
            raise BaselineError(
                f"fold {fold}: condition key {key!r} disagrees with its (drug, dose) pair {expected!r}"
            )

    counts_per_condition = np.bincount(inverse, minlength=n_conditions)
    if np.any(counts_per_condition == 0):
        raise BaselineError(f"fold {fold}: a training condition has no row")
    incidence = sp.csr_matrix(
        (
            np.ones(train_rows.size, dtype=np.float64),
            (np.arange(train_rows.size, dtype=np.int64), inverse),
        ),
        shape=(train_rows.size, n_conditions),
    )

    pool = control_pool or data_mod.CMControlPool(frozen)
    plate_matrix, plate_index = plate_control_matrix(pool, train_rows)

    n_genes = int(frozen.n_genes)
    condition_delta = np.empty((n_conditions, n_genes), dtype=np.float32)
    pooled_sum = np.zeros(n_genes, dtype=np.float64)
    for g0, g1, layer_block in iter_layer_blocks(
        frozen.cm_h5ad, train_rows, n_genes, block_genes=GENE_BLOCK
    ):
        control_block = plate_matrix[plate_index, g0:g1]
        delta_block = layer_block.astype(np.float64) - control_block.astype(np.float64)
        condition_block = np.asarray(incidence.T @ delta_block, dtype=np.float64)
        condition_block /= counts_per_condition[:, None]
        condition_delta[:, g0:g1] = condition_block.astype(np.float32)
        pooled_sum[g0:g1] = condition_block.sum(axis=0)
    pooled_mean = pooled_sum / float(n_conditions)

    plate = obs["plate"].astype(str).to_numpy()[train_rows]
    counts = {
        "fold": int(fold),
        "n_train_rows": int(train_rows.size),
        "n_train_conditions": int(unique.size),
        "n_train_drugs": int(np.unique(condition_drug).size),
        "n_train_plates": int(np.unique(plate).size),
        "train_dose_values": [float(d) for d in np.unique(condition_dose)],
        "n_conditions_per_drug": {
            str(drug_key): int(np.sum(condition_drug == drug_key))
            for drug_key in sorted(np.unique(condition_drug))
        },
        "memory_policy": {
            "layer_read": "dense log1p layer in gene blocks",
            "gene_block_genes": int(GENE_BLOCK),
            "condition_delta_dtype": "float32",
            "row_level_delta_materialised": False,
        },
    }
    return FitTable(
        fold=int(fold),
        n_genes=n_genes,
        train_rows=train_rows,
        row_conditions=inverse,
        unique_conditions=np.asarray(unique, dtype=object),
        condition_delta=condition_delta,
        condition_drug=np.asarray(condition_drug, dtype=object),
        condition_dose=np.asarray(condition_dose, dtype=np.float64),
        pooled_mean=np.asarray(pooled_mean, dtype=np.float64),
        counts=counts,
        train_rows_sha256=_sha256_rows(train_rows),
        train_delta_sha256=_sha256_array(condition_delta),
    )


# --------------------------------------------------------------------------- #
# test-side condition context
# --------------------------------------------------------------------------- #
@dataclass
class PredictionResult:
    delta: np.ndarray
    info: dict[str, Any]


def condition_context(
    frozen: data_mod.FrozenInput,
    source: data_mod.CMPredictSource,
    structures: data_mod.DrugStructures,
    rows: np.ndarray,
) -> dict[str, Any]:
    """Condition identities for *rows* in the evaluator's ``drug_key|dose`` order."""
    rows = np.asarray(rows, dtype=np.int64)
    keys = np.asarray(source.condition_key(rows), dtype=object)
    unique, inverse = np.unique(keys, return_inverse=True)
    inverse = np.asarray(inverse, dtype=np.int64)
    obs = frozen.cm_obs
    drug = obs["drug_key"].astype(str).to_numpy()[rows]
    dose = pd.to_numeric(obs["dose_um"], errors="coerce").to_numpy(dtype=np.float64)[rows]
    first = np.full(unique.size, -1, dtype=np.int64)
    for position, condition in enumerate(inverse):
        if first[condition] < 0:
            first[condition] = position
    if np.any(first < 0):
        raise BaselineError("a scored condition has no row")
    condition_drug = drug[first]
    condition_dose = dose[first]
    for i, key in enumerate(unique):
        expected = f"{condition_drug[i]}|{condition_dose[i]:g}"
        if str(key) != expected:
            raise BaselineError(
                f"condition key {key!r} disagrees with its (drug, dose) pair {expected!r}"
            )
    smiles = [structures.smiles_for(drug_key) or "" for drug_key in condition_drug]
    if any(not s for s in smiles):
        missing = sorted({str(d) for d, s in zip(condition_drug, smiles) if not s})[:5]
        raise BaselineError(f"no SMILES for scored condition drug(s) {missing}")
    return {
        "rows": rows,
        "condition_keys": np.asarray(unique, dtype=object),
        "row_conditions": inverse,
        "condition_drug": np.asarray(condition_drug, dtype=object),
        "condition_dose": np.asarray(condition_dose, dtype=np.float64),
        "condition_smiles": np.asarray(smiles, dtype=object),
        "inverse": inverse,
    }


def ridge_features(smiles: Sequence[str], doses: Sequence[float], n_bits: int) -> np.ndarray:
    """The study's dose encoding: 1024-bit rFCFP4 × log10(dose_µM + 1)."""
    smiles = np.asarray([str(s) for s in smiles], dtype=object)
    doses = np.asarray(doses, dtype=np.float64)
    if smiles.size != doses.size:
        raise BaselineError("smiles and doses must have the same length")
    unique, inverse = np.unique(smiles, return_inverse=True)
    fingerprints = data_mod.morgan_fingerprints(unique, int(n_bits)).astype(np.float64)
    scale = data_mod.dose_scale(doses).astype(np.float64)
    return fingerprints[np.asarray(inverse, dtype=np.int64)] * scale[:, None]


# --------------------------------------------------------------------------- #
# fits
# --------------------------------------------------------------------------- #
def fit_baseline(
    name: str,
    table: FitTable,
    structures: data_mod.DrugStructures,
    frozen: data_mod.FrozenInput,
    alpha: float = RIDGE_ALPHA,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Fit one baseline on the fold's training conditions; return (params, record)."""
    if name not in BASELINES:
        raise BaselineError(f"unknown baseline {name!r}; expected one of {list(BASELINES)}")

    if name == "no_change":
        params = {"delta": np.zeros(int(table.n_genes), dtype=np.float32)}
        record = {
            "baseline": name,
            "role": "primary",
            "rule": "Δ = 0 for every well; predicted expression = the fold's c*",
            "n_parameters": int(table.n_genes),
        }
    elif name == "delta_mean":
        params = {"pooled": np.asarray(table.pooled_mean, dtype=np.float32)}
        record = {
            "baseline": name,
            "role": "primary",
            "rule": (
                "within-condition replicate mean first, then equal-weight arithmetic mean of "
                "the fold's training condition means; a test condition whose exact "
                "drug_key|dose key also occurs in training uses that training condition's "
                "mean, otherwise the pooled training mean (documented fallback)"
            ),
            "pooling": "equal weight per training condition (drug_key|dose_µM), not per well",
            "lookup": "exact drug_key|dose_µM training-condition match, else pooled",
            "n_parameters": int(table.n_genes),
            "n_train_conditions": int(table.n_conditions),
        }
    elif name == "delta_mean_dose":
        dose_keys = np.asarray([f"{d:g}" for d in table.condition_dose], dtype=object)
        unique_doses = np.asarray(sorted(set(dose_keys.tolist())), dtype=object)
        means = np.stack(
            [
                table.condition_delta[dose_keys == dose].astype(np.float64).mean(axis=0)
                for dose in unique_doses
            ]
        ).astype(np.float32)
        params = {
            "doses": unique_doses.astype(str),
            "means": means,
            "pooled": np.asarray(table.pooled_mean, dtype=np.float32),
        }
        record = {
            "baseline": name,
            "role": "auxiliary",
            "rule": (
                "training condition means pooled per dose (:g µM); an unseen dose falls back "
                "to the pooled training mean"
            ),
            "doses": [str(d) for d in unique_doses],
            "n_conditions_per_dose": {
                str(d): int(np.sum(dose_keys == d)) for d in unique_doses
            },
            "n_parameters": int(len(unique_doses) * table.n_genes + table.n_genes),
        }
    else:  # ridge_fp
        from sklearn.linear_model import Ridge

        if structures is None:
            raise BaselineError("ridge_fp needs the frozen drug structures (SMILES)")
        n_bits = int(frozen.structure["drug_dim"])
        condition_smiles = [
            structures.smiles_for(drug_key) or "" for drug_key in table.condition_drug
        ]
        missing = sorted(
            {str(d) for d, s in zip(table.condition_drug, condition_smiles) if not s}
        )
        if missing:
            raise BaselineError(f"ridge_fp: no SMILES for training drug(s) {missing[:5]}")
        features = ridge_features(condition_smiles, table.condition_dose, n_bits)
        model = Ridge(alpha=float(alpha), fit_intercept=True, solver="cholesky")
        model.fit(features, table.condition_delta)
        params = {
            "coef": np.asarray(model.coef_, dtype=np.float32),
            "intercept": np.asarray(model.intercept_, dtype=np.float32),
        }
        record = {
            "baseline": name,
            "role": "auxiliary",
            "rule": (
                "L2 ridge from the study dose encoding (rFCFP4 1024-bit × log10(dose_µM+1)) "
                "to the training condition-mean Δ"
            ),
            "alpha": float(alpha),
            "alpha_policy": "fixed a priori at the 03.1 default; never tuned on any split",
            "fingerprint": "Morgan r=2, useFeatures=True (rFCFP4)",
            "n_bits": int(n_bits),
            "n_features": int(n_bits),
            "fit_intercept": True,
            "solver": "cholesky",
            "n_parameters": int(n_bits * table.n_genes + table.n_genes),
        }
    record["fold"] = int(table.fold)
    record["n_train_rows"] = int(table.n_rows)
    record["n_train_conditions"] = int(table.n_conditions)
    record["train_rows_sha256"] = table.train_rows_sha256
    record["train_delta_sha256"] = table.train_delta_sha256
    record["params_sha256"] = arrays_sha256(params)
    return params, record


# --------------------------------------------------------------------------- #
# predictions
# --------------------------------------------------------------------------- #
def condition_predictions(
    name: str,
    params: Mapping[str, np.ndarray],
    table: FitTable,
    context: Mapping[str, Any],
    frozen: data_mod.FrozenInput,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Predict one Δ per scored **condition** (``n_conditions × n_genes`` float32).

    Bounded by construction: the retained matrix is condition-level.  A row-level
    ``n_test_rows × n_genes`` array is only produced by the writer, one row block
    at a time.
    """
    n_genes = int(frozen.n_genes)
    n_conditions = int(context["condition_keys"].size)
    predicted = np.empty((n_conditions, n_genes), dtype=np.float32)

    if name == "no_change":
        predicted.fill(0.0)
        info: dict[str, Any] = {"rule": "zero Δ"}
    elif name == "delta_mean":
        pooled = np.asarray(params["pooled"], dtype=np.float32)
        if pooled.shape != (n_genes,):
            raise BaselineError(f"delta_mean pooled shape {pooled.shape} is not ({n_genes},)")
        predicted[:] = pooled[None, :]
        lookup = table.condition_lookup()
        matched = 0
        for i, key in enumerate(context["condition_keys"]):
            index = lookup.get(str(key))
            if index is not None:
                predicted[i] = table.condition_delta[index].astype(np.float32)
                matched += 1
        info = {
            "pooling": "training condition mean, else pooled training mean",
            "n_conditions_seen_in_train": int(matched),
            "fallback": "pooled training mean" if matched < n_conditions else "none",
        }
    elif name == "delta_mean_dose":
        pooled = np.asarray(params["pooled"], dtype=np.float32)
        means = np.asarray(params["means"], dtype=np.float32)
        doses = np.asarray(params["doses"]).astype(str)
        if means.shape != (doses.size, n_genes):
            raise BaselineError(
                f"delta_mean_dose means shape {means.shape} is not ({doses.size}, {n_genes})"
            )
        predicted[:] = pooled[None, :]
        by_dose = {str(dose): i for i, dose in enumerate(doses)}
        missing = 0
        for i, dose in enumerate(context["condition_dose"]):
            index = by_dose.get(f"{float(dose):g}")
            if index is None:
                missing += 1
            else:
                predicted[i] = means[index]
        info = {
            "pooling": "per training dose",
            "n_conditions_with_unseen_dose": int(missing),
            "fallback": "pooled training mean",
        }
    elif name == "ridge_fp":
        coef = np.asarray(params["coef"], dtype=np.float32)
        if coef.ndim == 1:  # single-gene fit; sklearn collapses the target axis
            coef = coef[None, :]
        intercept = np.asarray(params["intercept"], dtype=np.float32).reshape(-1)
        n_bits = int(coef.shape[1])
        if coef.shape[0] != n_genes or intercept.size != n_genes:
            raise BaselineError(
                f"ridge_fp coefficient shape {coef.shape} / intercept {intercept.shape} "
                f"does not match ({n_genes}, {n_bits})"
            )
        features = ridge_features(
            context["condition_smiles"], context["condition_dose"], n_bits
        ).astype(np.float64)
        for g0 in range(0, n_genes, GENE_BLOCK):
            g1 = min(g0 + GENE_BLOCK, n_genes)
            block = features @ coef[g0:g1].astype(np.float64).T
            block += intercept[g0:g1].astype(np.float64)[None, :]
            predicted[:, g0:g1] = block.astype(np.float32)
        info = {"model": "ridge", "n_features": n_bits}
    else:  # pragma: no cover - fit validated the name
        raise BaselineError(f"unknown baseline {name!r}")

    if predicted.shape != (n_conditions, n_genes):
        raise BaselineError(f"baseline {name} produced Δ shape {predicted.shape}")
    if not np.all(np.isfinite(predicted)):
        raise BaselineError(f"baseline {name} produced non-finite Δ")
    return predicted, info


def predict_baseline(
    name: str,
    params: Mapping[str, np.ndarray],
    table: FitTable,
    context: Mapping[str, Any],
    frozen: data_mod.FrozenInput,
) -> PredictionResult:
    """Row-level Δ for the scored rows.

    Test/small-input helper only: it materialises ``n_rows × n_genes``.  The
    five-fold run path uses :func:`condition_predictions` + the blocked writer so
    that array is never built in memory.
    """
    inverse = np.asarray(context["row_conditions"], dtype=np.int64)
    predicted, info = condition_predictions(name, params, table, context, frozen)
    delta = np.ascontiguousarray(predicted[inverse], dtype=np.float32)
    if delta.shape != (inverse.size, int(frozen.n_genes)):
        raise BaselineError(f"baseline {name} produced Δ shape {delta.shape}")
    return PredictionResult(delta=delta, info=info)


# --------------------------------------------------------------------------- #
# prediction directories in the pipeline format
# --------------------------------------------------------------------------- #
def write_prediction_dir(
    out_dir: Path,
    *,
    baseline: str,
    role: str,
    fold: int,
    meta: pd.DataFrame,
    c_star: np.ndarray,
    fit_record: Mapping[str, Any],
    frozen: data_mod.FrozenInput,
    command: Sequence[str] = (),
    repeats: int = 2,
    delta: np.ndarray | None = None,
    condition_delta: np.ndarray | None = None,
    row_conditions: np.ndarray | None = None,
) -> dict[str, Any]:
    """Write ``pred_expr.npy`` / ``delta.npy`` / ``meta.csv`` / ``manifest.json``.

    ``pred_expr = c* + Δ`` and the stored Δ is recomputed as ``pred_expr − c*``
    so the evaluator's baseline check passes exactly; the caller's Δ is also
    recomputed and required to be byte-identical (the deterministic-prediction
    policy, mirroring ``ptlib.predict``).

    Memory (ticket #8): pass ``condition_delta`` + ``row_conditions`` from
    :func:`condition_predictions` and the two ``.npy`` files are created as
    ``np.lib.format`` memmaps and filled in ``ROW_BLOCK`` row blocks — no
    ``n_rows × n_genes`` array is held in memory.  ``delta=`` (a row-level array)
    is kept for the evaluator-sanity tests, which are fixture-sized.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n_rows = int(len(meta))
    n_genes = int(frozen.n_genes)
    if n_rows == 0:
        raise BaselineError("cannot write an empty prediction directory")
    c_star64 = np.ascontiguousarray(np.asarray(c_star, dtype=np.float32)).astype(np.float64)
    if c_star64.shape != (n_genes,):
        raise BaselineError(f"c* shape {c_star64.shape} is not ({n_genes},)")
    if (delta is None) == (condition_delta is None):
        raise BaselineError(
            "pass exactly one of delta (row-level) or condition_delta (+ row_conditions)"
        )

    if condition_delta is not None:
        condition_delta = np.asarray(condition_delta, dtype=np.float32)
        if condition_delta.ndim != 2 or condition_delta.shape[1] != n_genes:
            raise BaselineError(
                f"condition Δ shape {condition_delta.shape} is not (*, {n_genes})"
            )
        row_conditions = np.asarray(row_conditions, dtype=np.int64)
        if row_conditions.shape != (n_rows,):
            raise BaselineError(
                f"row_conditions shape {row_conditions.shape} is not ({n_rows},)"
            )
        if row_conditions.size and (
            row_conditions.min() < 0 or row_conditions.max() >= condition_delta.shape[0]
        ):
            raise BaselineError("row_conditions index outside the condition Δ matrix")

        def delta_block(start: int, stop: int) -> np.ndarray:
            return condition_delta[row_conditions[start:stop]]

    else:
        delta = np.asarray(delta, dtype=np.float32)
        if delta.shape != (n_rows, n_genes):
            raise BaselineError(
                f"prediction Δ shape {delta.shape} does not match meta rows {n_rows} "
                f"and {n_genes} genes"
            )

        def delta_block(start: int, stop: int) -> np.ndarray:
            return delta[start:stop]

    pred_path = out_dir / "pred_expr.npy"
    stored_path = out_dir / "delta.npy"
    predicted_expr = np.lib.format.open_memmap(
        pred_path, mode="w+", dtype=np.float32, shape=(n_rows, n_genes)
    )
    stored_delta = np.lib.format.open_memmap(
        stored_path, mode="w+", dtype=np.float32, shape=(n_rows, n_genes)
    )
    baseline_roundtrip = 0.0
    identical = True
    extra_repeats = max(1, int(repeats)) - 1
    try:
        for start in range(0, n_rows, ROW_BLOCK):
            stop = min(start + ROW_BLOCK, n_rows)
            block = np.ascontiguousarray(delta_block(start, stop), dtype=np.float32)
            if block.shape != (stop - start, n_genes):
                raise BaselineError(
                    f"prediction Δ block {block.shape} does not match ({stop - start}, {n_genes})"
                )
            block64 = block.astype(np.float64)
            pred_block = (c_star64[None, :] + block64).astype(np.float32)
            stored_block = (pred_block.astype(np.float64) - c_star64).astype(np.float32)
            predicted_expr[start:stop] = pred_block
            stored_delta[start:stop] = stored_block
            roundtrip = float(
                np.max(
                    np.abs(
                        pred_block.astype(np.float64)
                        - stored_block.astype(np.float64)
                        - c_star64
                    )
                )
            )
            baseline_roundtrip = max(baseline_roundtrip, roundtrip)
            for _ in range(extra_repeats):
                again = (c_star64[None, :] + block64).astype(np.float32)
                if again.tobytes() != pred_block.tobytes():
                    identical = False
    finally:
        predicted_expr.flush()
        stored_delta.flush()
        del predicted_expr, stored_delta

    if baseline_roundtrip > 1e-5:
        raise BaselineError(
            f"prediction baseline round-trip differs from c* by {baseline_roundtrip:.3g}"
        )
    if not identical:
        raise BaselineError(
            f"baseline {baseline}: repeated prediction is not byte-identical "
            "(deterministic prediction policy violated)"
        )
    meta.to_csv(out_dir / "meta.csv", index=False)

    rows = np.asarray(meta["row"], dtype=np.int64)
    pred_expr_sha256 = manifest_mod.sha256_file(pred_path)
    delta_sha256 = manifest_mod.sha256_file(stored_path)
    meta_sha256 = manifest_mod.sha256_file(out_dir / "meta.csv")
    payload = {
        "schema": "stage2/06.8_baselines/prediction/v1",
        "kind": "baselines_prediction",
        "baseline": baseline,
        "role": role,
        "stage": "finetune",
        "scored_role": "test",
        "fold": int(fold),
        "rows": int(n_rows),
        "n_genes": n_genes,
        "row_ids_sha256": _sha256_rows(rows),
        "c_star_sha256": _sha256_array(np.asarray(c_star, dtype=np.float32)),
        "predictions": {
            # File-content sha256 (the evaluator hashes the same files); the
            # arrays are written through memmaps and never hashed in memory.
            "pred_expr_sha256": pred_expr_sha256,
            "delta_sha256": delta_sha256,
            "meta_sha256": meta_sha256,
            "pred_expr_shape": [int(n_rows), int(n_genes)],
            "input_baseline_max_abs_diff": baseline_roundtrip,
        },
        "determinism": {
            "noise": "none",
            "repeats": int(max(1, int(repeats))),
            "repeated_inference_identical": bool(identical),
            "seed": None,
            "note": "baselines are closed-form/deterministic; there is no model seed",
        },
        "config": {
            "baseline": baseline,
            "role": role,
            "ridge_alpha": RIDGE_ALPHA if baseline == "ridge_fp" else None,
            "delta": "prediction − the fold's fixed c*",
            "truth_delta": "treated log1p_cp10k − per-plate matched DMSO mean (evaluator only)",
            "replicate_aggregation": "within condition first (evaluator)",
            "memory": {
                "row_block_rows": int(ROW_BLOCK),
                "gene_block_genes": int(GENE_BLOCK),
                "prediction_arrays": (
                    "np.lib.format.open_memmap written in row blocks; no dense "
                    "n_rows × n_genes array in memory"
                ),
            },
        },
        "fit": dict(fit_record),
        "data": frozen.record(),
        "git": manifest_mod.git_info(ROOT),
        "command": list(command),
        "created_utc": manifest_mod.now_utc(),
    }
    # The shared evaluator reads this key; baselines have no checkpoint.
    payload["checkpoint"] = None
    manifest_mod.write_json(out_dir / "manifest.json", payload)
    return {
        "dir": str(out_dir),
        "rows": int(n_rows),
        "row_ids_sha256": payload["row_ids_sha256"],
        "c_star_sha256": payload["c_star_sha256"],
        "prediction_sha256": pred_expr_sha256,
        "delta_sha256": delta_sha256,
        "meta_sha256": meta_sha256,
        "input_baseline_max_abs_diff": baseline_roundtrip,
        "repeated_inference_identical": bool(identical),
    }


def score_prediction_dir(
    frozen: data_mod.FrozenInput,
    pred_dir: Path,
    eval_dir: Path,
    fold: int,
    *,
    hash_cache: Path | str | None = None,
    table_cache: Path | str | None = None,
    command: Sequence[str] = (),
) -> dict[str, Any]:
    """Score a baseline prediction directory with the shared evaluator."""
    options = evaluate_mod.EvaluateOptions(
        stage="finetune",
        pred_dir=Path(pred_dir),
        out_dir=Path(eval_dir),
        data_dir=frozen.data_dir,
        role="test",
        fold=int(fold),
        hash_cache=Path(hash_cache) if hash_cache is not None else None,
        table_cache=Path(table_cache) if table_cache is not None else None,
        command=list(command),
    )
    return evaluate_mod.run_evaluation(options)


# --------------------------------------------------------------------------- #
# metrics bookkeeping
# --------------------------------------------------------------------------- #
def headline_metrics(aggregate: Mapping[str, Any], k: int) -> dict[str, Any]:
    """Extract the headline (k=200) macro metrics from an evaluator aggregate."""
    by_k = aggregate["by_k"][str(int(k))]
    regression = aggregate["regression"]

    def _stat(block: Mapping[str, Any]) -> dict[str, float]:
        return {
            "mean": float(block["mean"]),
            "sd": float(block["sd"]) if block["sd"] is not None else math.nan,
            "nan": int(block["nan"]),
        }

    macro = _stat(by_k["macro"])
    up = _stat(by_k["up"])
    down = _stat(by_k["down"])
    rmse = _stat(regression["rmse"])
    pcc = _stat(regression["pcc"])
    spearman = _stat(regression["spearman"])
    return {
        "k": int(k),
        "auprc": macro["mean"],
        "auprc_sd_conditions": macro["sd"],
        "auprc_nan": macro["nan"],
        "auprc_up": up["mean"],
        "auprc_down": down["mean"],
        "rmse": rmse["mean"],
        "pcc": pcc["mean"],
        "spearman": spearman["mean"],
        "n_conditions": int(aggregate["n_conditions"]),
        "n_rows": int(aggregate["n_rows"]),
        "nan_rmse": int(rmse["nan"]),
        "nan_pcc": int(pcc["nan"]),
        "nan_spearman": int(spearman["nan"]),
    }


def summarise_folds(per_fold: Sequence[Mapping[str, Any]], baselines: Sequence[str]) -> list[dict[str, Any]]:
    """Row records for ``summary_table.csv``: per fold plus macro-over-folds."""
    metrics = ("auprc", "auprc_up", "auprc_down", "rmse", "pcc", "spearman")
    records: list[dict[str, Any]] = []
    for baseline in baselines:
        rows = [r for r in per_fold if r["baseline"] == baseline]
        if not rows:
            continue
        for row in sorted(rows, key=lambda item: item["fold"]):
            metrics_row = row["metrics"]
            record: dict[str, Any] = {
                "scope": "per_fold",
                "fold": int(row["fold"]),
                "baseline": baseline,
                "role": row["role"],
                "n_test_rows": int(metrics_row["n_rows"]),
                "n_test_conditions": int(metrics_row["n_conditions"]),
            }
            for metric in metrics:
                record[metric] = float(metrics_row[metric])
                record[f"{metric}_sd"] = math.nan
            record["auprc_nan_conditions"] = int(metrics_row["auprc_nan"])
            records.append(record)

        record = {
            "scope": "macro_over_folds",
            "fold": "ALL",
            "baseline": baseline,
            "role": rows[0]["role"],
            "n_test_rows": int(sum(r["metrics"]["n_rows"] for r in rows)),
            "n_test_conditions": int(sum(r["metrics"]["n_conditions"] for r in rows)),
        }
        for metric in metrics:
            values = np.asarray([float(r["metrics"][metric]) for r in rows], dtype=np.float64)
            finite = values[np.isfinite(values)]
            record[metric] = float(finite.mean()) if finite.size else math.nan
            record[f"{metric}_sd"] = (
                float(np.std(finite, ddof=1)) if finite.size > 1 else math.nan
            )
        record["auprc_nan_conditions"] = int(sum(r["metrics"]["auprc_nan"] for r in rows))
        records.append(record)
    return records


# --------------------------------------------------------------------------- #
# one fold
# --------------------------------------------------------------------------- #
def run_fold(
    frozen: data_mod.FrozenInput,
    fold_table: data_mod.FoldTable,
    fold: int,
    out_dir: Path,
    baselines: Sequence[str],
    *,
    alpha: float = RIDGE_ALPHA,
    hash_cache: Path | str | None = None,
    table_cache: Path | str | None = None,
    force: bool = False,
    command: Sequence[str] = (),
) -> dict[str, Any]:
    """Fit every requested baseline on *fold*'s training responses and score it."""
    started = time.perf_counter()
    fold_dir = Path(out_dir) / "folds" / f"fold{int(fold)}"
    fold_dir.mkdir(parents=True, exist_ok=True)
    structures = frozen.structures()
    source = data_mod.CMPredictSource(frozen, fold_table, int(fold), structures)
    control_pool = data_mod.CMControlPool(frozen)
    table = build_fit_table(
        frozen, fold_table, int(fold), structures=structures,
        control_pool=control_pool, source=source,
    )
    test_rows = np.asarray(source.rows_for("test"), dtype=np.int64)
    if test_rows.size == 0:
        raise BaselineError(f"fold {fold}: no test treated rows")
    valid_rows = np.asarray(source.rows_for("valid"), dtype=np.int64)
    obs = frozen.cm_obs
    test_drug = obs["drug_key"].astype(str).to_numpy()[test_rows]
    meta = source.meta(test_rows)
    c_star = np.asarray(source.c_star, dtype=np.float32)
    context = condition_context(frozen, source, structures, test_rows)
    k = int(frozen.config["metric"]["k"])

    records: list[dict[str, Any]] = []
    for name in baselines:
        role = "primary" if name in PRIMARY_BASELINES else "auxiliary"
        fit_dir = fold_dir / "fits" / name
        pred_dir = fold_dir / name / "pred"
        eval_dir = fold_dir / name / "eval"
        aggregate_file = eval_dir / "aggregate.json"
        if (
            not force
            and aggregate_file.exists()
            and (pred_dir / "manifest.json").exists()
        ):
            aggregate = manifest_mod.read_json(aggregate_file)
            pred_manifest = manifest_mod.read_json(pred_dir / "manifest.json")
            fit_record = dict(pred_manifest.get("fit") or {})
            prediction = {
                "dir": str(pred_dir),
                "rows": int(pred_manifest.get("rows", len(meta))),
                "row_ids_sha256": pred_manifest.get("row_ids_sha256"),
                "c_star_sha256": pred_manifest.get("c_star_sha256"),
                "prediction_sha256": (pred_manifest.get("predictions") or {}).get(
                    "pred_expr_sha256"
                ),
                "delta_sha256": (pred_manifest.get("predictions") or {}).get("delta_sha256"),
                "meta_sha256": (pred_manifest.get("predictions") or {}).get("meta_sha256"),
                "reused": True,
            }
        else:
            params, fit_record = fit_baseline(
                name, table, structures, frozen, alpha=alpha
            )
            manifest_mod.write_json(fit_dir / "fit_record.json", fit_record)
            np.savez_compressed(fit_dir / "params.npz", **params)
            condition_delta, prediction_info = condition_predictions(
                name, params, table, context, frozen
            )
            fit_record = {**fit_record, "prediction_rule": prediction_info}
            written = write_prediction_dir(
                pred_dir,
                baseline=name,
                role=role,
                fold=int(fold),
                meta=meta,
                c_star=c_star,
                condition_delta=condition_delta,
                row_conditions=context["row_conditions"],
                fit_record=fit_record,
                frozen=frozen,
                command=list(command),
            )
            del condition_delta
            result = score_prediction_dir(
                frozen,
                pred_dir,
                eval_dir,
                int(fold),
                hash_cache=hash_cache,
                table_cache=table_cache,
                command=list(command),
            )
            aggregate = result["aggregate"]
            prediction = {**written, "reused": False, "rule": prediction_info}

        metrics = headline_metrics(aggregate, k)
        eval_manifest = manifest_mod.read_json(eval_dir / "manifest.json")
        records.append(
            {
                "fold": int(fold),
                "baseline": name,
                "role": role,
                "metrics": metrics,
                "prediction": prediction,
                "fit": fit_record,
                "eval_dir": str(eval_dir),
                "pred_dir": str(pred_dir),
                "eval_manifest_sha256": manifest_mod.sha256_file(eval_dir / "manifest.json"),
                "truth": eval_manifest.get("truth"),
            }
        )

    per_condition = pd.read_csv(fold_dir / baselines[0] / "eval" / "per_condition.csv")
    coverage = {
        "fold": int(fold),
        "n_train_rows": int(table.n_rows),
        "n_train_conditions": int(table.n_conditions),
        "n_train_drugs": int(table.counts["n_train_drugs"]),
        "n_valid_rows": int(valid_rows.size),
        "n_test_rows": int(test_rows.size),
        "n_test_conditions": int(context["condition_keys"].size),
        "n_test_drugs": int(np.unique(test_drug).size),
        "n_test_conditions_seen_in_train": int(
            sum(
                str(key) in table.condition_lookup()
                for key in context["condition_keys"]
            )
        ),
        "n_plates_train": int(table.counts["n_train_plates"]),
        "train_row_sha256": table.train_rows_sha256,
        "test_row_sha256": _sha256_rows(test_rows),
        "c_star_sha256": _sha256_array(c_star),
        "n_scored_conditions": int(len(per_condition)),
        "n_scored_replicates": int(per_condition["n_replicates"].sum()),
        "scored_conditions_unique": bool(not per_condition["condition"].duplicated().any()),
        "gene_axis_sha256": _sha256_array(np.asarray(frozen.var_names, dtype=str)),
        "test_rows_disjoint_from_train": bool(
            np.intersect1d(test_rows, table.train_rows).size == 0
        ),
    }
    fold_record = {
        "schema": SCHEMA,
        "kind": "baselines_fold",
        "fold": int(fold),
        "k": int(k),
        "coverage": coverage,
        "fit_counts": table.counts,
        "baselines": records,
        "git": manifest_mod.git_info(ROOT),
        "timings": {"wall_s": time.perf_counter() - started},
        "peak_rss_mb": manifest_mod.peak_rss_mb(),
    }
    manifest_mod.write_json(fold_dir / "fold_record.json", fold_record)
    return {"coverage": coverage, "baselines": records, "fold_record": fold_record}


# --------------------------------------------------------------------------- #
# all folds + the tracked record
# --------------------------------------------------------------------------- #
def _coverage_frame(folds: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    columns = [
        "fold", "n_train_rows", "n_train_conditions", "n_train_drugs", "n_valid_rows",
        "n_test_rows", "n_test_conditions", "n_test_drugs",
        "n_test_conditions_seen_in_train", "n_plates_train", "n_scored_conditions",
        "n_scored_replicates", "scored_conditions_unique", "test_rows_disjoint_from_train",
        "train_row_sha256", "test_row_sha256", "c_star_sha256", "gene_axis_sha256",
    ]
    return pd.DataFrame([{column: fold["coverage"][column] for column in columns} for fold in folds])


def _per_drug_dose_frame(folds: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    frames = []
    for fold in folds:
        for record in fold["baselines"]:
            eval_dir = Path(record["eval_dir"])
            frame = pd.read_csv(eval_dir / "per_condition.csv")
            # The evaluator persists the true-set size only for "up"; the down
            # count is recoverable from the same deltas.npz it wrote (same
            # condition order).  #9 needs both for the direction split.
            payload = np.load(eval_dir / "deltas.npz")
            truth = np.asarray(payload["true_delta"], dtype=np.float64)
            k = int(record["metrics"]["k"])
            n_pos = (truth > 0).sum(axis=1)
            n_neg = (truth < 0).sum(axis=1)
            expected_up = np.minimum(n_pos, k)
            if not np.array_equal(expected_up, frame[f"n_up_{k}"].to_numpy(dtype=np.int64)):
                raise BaselineError(
                    f"per_condition.csv up-set sizes disagree with deltas.npz in {eval_dir}"
                )
            frame[f"n_down_{k}"] = np.minimum(n_neg, k).astype(np.int64)
            frame.insert(0, "fold", int(fold["coverage"]["fold"]))
            frame.insert(1, "baseline", record["baseline"])
            frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    condition = combined["condition"].astype(str)
    combined.insert(
        2, "drug_key", condition.str.rsplit("|", n=1).str[0]
    )
    combined.insert(
        3, "dose_um", pd.to_numeric(condition.str.rsplit("|", n=1).str[1], errors="coerce")
    )
    columns = [
        "fold", "baseline", "drug_key", "dose_um", "condition", "n_replicates",
        "auprc_200", "auprc_up_200", "auprc_down_200", "n_up_200", "n_down_200",
        "rmse", "pcc", "spearman", "pred_constant", "true_constant",
        "regression_nan_reason",
    ]
    missing = [column for column in columns if column not in combined.columns]
    if missing:
        raise BaselineError(f"per_condition.csv is missing {missing}")
    return combined[columns]


def _metrics_frame(folds: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    frames = []
    for fold in folds:
        for record in fold["baselines"]:
            frame = pd.read_csv(Path(record["eval_dir"]) / "metrics.csv")
            frame.insert(0, "fold", int(fold["coverage"]["fold"]))
            frame.insert(1, "baseline", record["baseline"])
            frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _check_fold_invariants(
    frozen: data_mod.FrozenInput,
    fold_table: data_mod.FoldTable,
    folds: Sequence[Mapping[str, Any]],
    baselines: Sequence[str],
) -> dict[str, Any]:
    n_folds = int(frozen.config["n_folds"])
    coverage = [fold["coverage"] for fold in folds]
    main_drugs = set(fold_table.drugs.tolist())
    test_drug_occurrences: dict[str, list[int]] = {}
    for fold in range(n_folds):
        for drug in fold_table.drugs_for(fold, "test"):
            test_drug_occurrences.setdefault(str(drug), []).append(fold)
    invariants = {
        "n_folds": n_folds,
        "n_main_drugs": int(len(main_drugs)),
        "every_main_drug_tested_exactly_once": bool(
            len(test_drug_occurrences) == len(main_drugs)
            and all(len(v) == 1 for v in test_drug_occurrences.values())
        ),
        "test_rows_disjoint_from_train": bool(
            all(row["test_rows_disjoint_from_train"] for row in coverage)
        ),
        "no_test_condition_seen_in_train": bool(
            all(row["n_test_conditions_seen_in_train"] == 0 for row in coverage)
        ),
        "every_test_condition_scored_once": bool(
            all(
                row["scored_conditions_unique"]
                and row["n_scored_conditions"] == row["n_test_conditions"]
                and row["n_scored_replicates"] == row["n_test_rows"]
                for row in coverage
            )
        ),
        "same_test_population_across_baselines": bool(
            all(
                len(
                    {
                        (
                            record["prediction"].get("row_ids_sha256"),
                            record["prediction"].get("meta_sha256"),
                            record["metrics"]["n_rows"],
                            record["metrics"]["n_conditions"],
                        )
                        for record in fold["baselines"]
                        if record["baseline"] in baselines
                    }
                )
                == 1
                for fold in folds
            )
        ),
        "c_star_identical_across_baselines": bool(
            all(
                len(
                    {
                        record["prediction"].get("c_star_sha256")
                        for record in fold["baselines"]
                        if record["baseline"] in baselines
                    }
                )
                == 1
                for fold in folds
            )
        ),
    }
    return invariants


def run_all(
    *,
    data_dir: Path | str,
    out_dir: Path | str,
    folds: Sequence[int] | None = None,
    baselines: Sequence[str] = BASELINES,
    alpha: float = RIDGE_ALPHA,
    hash_cache: Path | str | None = None,
    table_cache: Path | str | None = None,
    force: bool = False,
    command: Sequence[str] = (),
) -> dict[str, Any]:
    """Run every fold, score the baselines, and write the tracked record."""
    started = time.perf_counter()
    started_utc = manifest_mod.now_utc()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frozen = data_mod.open_frozen_input(
        data_dir, hash_cache=hash_cache, table_cache=table_cache
    )
    fold_table = data_mod.load_folds(frozen)
    fold_numbers = (
        list(range(int(frozen.config["n_folds"]))) if folds is None else [int(f) for f in folds]
    )
    for fold in fold_numbers:
        if fold < 0 or fold >= int(frozen.config["n_folds"]):
            raise BaselineError(
                f"fold {fold} outside the frozen n_folds={frozen.config['n_folds']}"
            )
    used = [name for name in baselines]
    unknown = [name for name in used if name not in BASELINES]
    if unknown:
        raise BaselineError(f"unknown baseline(s) {unknown}; expected {list(BASELINES)}")

    fold_results = []
    for fold in fold_numbers:
        result = run_fold(
            frozen,
            fold_table,
            fold,
            out_dir,
            used,
            alpha=alpha,
            hash_cache=hash_cache,
            table_cache=table_cache,
            force=force,
            command=command,
        )
        fold_results.append(result)
        print(
            f"[06.8] fold {fold}: {result['coverage']['n_test_rows']} test rows, "
            f"{result['coverage']['n_test_conditions']} test conditions",
            flush=True,
        )

    summary = summarise_folds(
        [record for fold in fold_results for record in fold["baselines"]], used
    )
    summary_frame = pd.DataFrame(summary)
    metrics_frame = _metrics_frame(fold_results)
    per_drug_dose = _per_drug_dose_frame(fold_results)
    coverage_frame = _coverage_frame(fold_results)

    summary_frame.to_csv(out_dir / "summary_table.csv", index=False, float_format="%.10g")
    metrics_frame.to_csv(out_dir / "metrics.csv", index=False, float_format="%.10g")
    per_drug_dose.to_csv(out_dir / "per_drug_dose.csv", index=False, float_format="%.10g")
    coverage_frame.to_csv(out_dir / "fold_coverage.csv", index=False, float_format="%.10g")

    invariants = _check_fold_invariants(frozen, fold_table, fold_results, used)
    if not all(
        value
        for key, value in invariants.items()
        if isinstance(value, bool) and key != "no_test_condition_seen_in_train"
    ):
        raise BaselineError(f"fold-coverage invariants failed: {invariants}")
    if not invariants["no_test_condition_seen_in_train"]:
        print(
            "[06.8] warning: a test condition also appears in training; the documented "
            "per-condition lookup rule was applied (recorded in the per-fold records)",
            file=sys.stderr,
        )

    payload = {
        "schema": SCHEMA,
        "kind": "baselines",
        "ticket": 8,
        "status": "completed",
        "created_utc": started_utc,
        "finished_utc": manifest_mod.now_utc(),
        "pipeline_version": PIPELINE_VERSION,
        "git": manifest_mod.git_info(ROOT),
        "command": list(command),
        "data": frozen.record(),
        "input": {
            "dir": str(frozen.data_dir),
            "artifact_sha256": dict(frozen.artifact_hashes),
            "manifest_sha256": frozen.manifest_sha256,
            "locked_config": {
                "path": str(frozen.data_dir / "locked_config.json"),
                "sha256": manifest_mod.sha256_file(frozen.data_dir / "locked_config.json"),
                "mapping": frozen.config_notes.get("mapping", []),
                "defaulted": frozen.config_notes.get("defaulted", []),
            },
            "gene_axis": {
                "n_genes": int(frozen.n_genes),
                "sha256": _sha256_array(np.asarray(frozen.var_names, dtype=str)),
            },
        },
        "config": {
            "baselines": list(used),
            "primary_baselines": [name for name in used if name in PRIMARY_BASELINES],
            "auxiliary_baselines": [name for name in used if name in AUXILIARY_BASELINES],
            "ridge_alpha": float(alpha),
            "ridge_alpha_policy": "fixed a priori (03.1 default); never tuned",
            "dose_encoding": dict(frozen.config["dose_encoding"]),
            "folds": fold_numbers,
            "metric": {
                "name": "direction_auprc",
                "headline_k": int(frozen.config["metric"]["k"]),
                "k_list": [int(k) for k in metrics_mod.K_LIST if int(k) <= frozen.n_genes],
                "aggregation": (
                    "within-condition replicate mean first (condition = drug_key|dose_µM), "
                    "then per-direction AP, macro over directions, macro over conditions"
                ),
                "truth_delta": "treated log1p_cp10k − per-plate matched DMSO mean (evaluator only)",
                "auxiliary": ["Δ RMSE", "Pearson", "Spearman"],
            },
            "fit_rule": (
                "only the fold's role=train treated responses are fitted; validation/test "
                "treated rows are never read"
            ),
            "pooling": {
                "delta_mean": "equal weight per training condition (drug_key|dose_µM)",
                "delta_mean_dose": "per training dose",
                "lookup": "exact training-condition match, else painted training-pool fallback",
            },
            "prediction_dir_format": (
                "ptlib.predict format: pred_expr.npy, delta.npy, meta.csv (row ids), "
                "manifest.json; Δ = prediction − c*"
            ),
            "memory": {
                "gene_block_genes": int(GENE_BLOCK),
                "row_block_rows": int(ROW_BLOCK),
                "row_level_delta_materialised": False,
                "prediction_arrays": "np.lib.format memmaps filled in row blocks",
                "n_jobs": "single process; no parallel fits",
                "policy": (
                    "no dense observation × n_genes array: the dense log1p layer is read "
                    "and aggregated in gene blocks, and predictions are written in row "
                    "blocks (ticket #8 memory rule)"
                ),
            },
        },
        "folds": {
            str(result["coverage"]["fold"]): result["fold_record"]["coverage"]
            for result in fold_results
        },
        "invariants": invariants,
        "outputs": {
            "root": str(out_dir),
            "summary_table": "summary_table.csv",
            "metrics": "metrics.csv",
            "per_drug_dose": "per_drug_dose.csv",
            "fold_coverage": "fold_coverage.csv",
            "baseline_runs": "folds/fold<k>/<baseline>/{pred,eval}/",
            "fit_params": "folds/fold<k>/fits/<baseline>/params.npz (local, not tracked)",
        },
        "record_sha256": {
            "summary_table": manifest_mod.sha256_file(out_dir / "summary_table.csv"),
            "metrics": manifest_mod.sha256_file(out_dir / "metrics.csv"),
            "per_drug_dose": manifest_mod.sha256_file(out_dir / "per_drug_dose.csv"),
            "fold_coverage": manifest_mod.sha256_file(out_dir / "fold_coverage.csv"),
        },
        "environment": _environment(),
        "timings": {
            "wall_s": time.perf_counter() - started,
            "per_fold_s": [
                result["fold_record"]["timings"]["wall_s"] for result in fold_results
            ],
        },
        "peak_rss_mb": manifest_mod.peak_rss_mb(),
        "warnings": (
            []
            if invariants["no_test_condition_seen_in_train"]
            else ["a test condition occurred in training; per-condition lookup was applied"]
        ),
    }
    manifest_mod.write_json(out_dir / "manifest.json", payload)
    print(
        json.dumps(
            {
                "out_dir": str(out_dir),
                "folds": fold_numbers,
                "baselines": used,
                "summary_rows": int(len(summary_frame)),
                "per_drug_dose_rows": int(len(per_drug_dose)),
                "n_test_rows_total": int(
                    summary_frame[summary_frame["scope"] == "per_fold"]["n_test_rows"].sum()
                ),
                "invariants": invariants,
            },
            indent=2,
            default=str,
        )
    )
    return {"manifest": payload, "summary": summary_frame, "out_dir": str(out_dir)}


# --------------------------------------------------------------------------- #
# contract check (Phase A / frozen-input readiness)
# --------------------------------------------------------------------------- #
def check_input(
    *,
    data_dir: Path | str,
    folds: Sequence[int] | None = None,
    hash_cache: Path | str | None = None,
    table_cache: Path | str | None = None,
) -> dict[str, Any]:
    """Verify the frozen input is readable and the baseline contract holds.

    Reads only metadata + the training treated rows (to count condition-key
    collisions).  Writes nothing.
    """
    frozen = data_mod.open_frozen_input(
        data_dir, hash_cache=hash_cache, table_cache=table_cache
    )
    fold_table = data_mod.load_folds(frozen)
    fold_numbers = (
        list(range(int(frozen.config["n_folds"]))) if folds is None else [int(f) for f in folds]
    )
    structures = frozen.structures()
    control_pool = data_mod.CMControlPool(frozen)
    per_fold = {}
    for fold in fold_numbers:
        source = data_mod.CMPredictSource(frozen, fold_table, int(fold), structures)
        table = build_fit_table(
            frozen, fold_table, int(fold), structures=structures,
            control_pool=control_pool, source=source,
        )
        test_rows = np.asarray(source.rows_for("test"), dtype=np.int64)
        context = condition_context(frozen, source, structures, test_rows)
        collisions = int(
            sum(str(key) in table.condition_lookup() for key in context["condition_keys"])
        )
        per_fold[str(fold)] = {
            **table.counts,
            "n_valid_rows": int(np.asarray(source.rows_for("valid")).size),
            "n_test_rows": int(test_rows.size),
            "n_test_conditions": int(context["condition_keys"].size),
            "n_test_drugs": int(
                np.unique(frozen.cm_obs["drug_key"].astype(str).to_numpy()[test_rows]).size
            ),
            "n_test_conditions_seen_in_train": collisions,
            "c_star_sha256": _sha256_array(np.asarray(source.c_star, dtype=np.float32)),
        }
    report = {
        "status": "ok",
        "data_dir": str(frozen.data_dir),
        "n_genes": int(frozen.n_genes),
        "n_folds": int(frozen.config["n_folds"]),
        "headline_k": int(frozen.config["metric"]["k"]),
        "artifact_hashes_verified": len(frozen.hash_record.get("verified", {}) or {}),
        "folds": per_fold,
        "baselines": list(BASELINES),
        "primary_baselines": list(PRIMARY_BASELINES),
        "auxiliary_baselines": list(AUXILIARY_BASELINES),
        "ridge_alpha": float(RIDGE_ALPHA),
        "memory": {
            "gene_block_genes": int(GENE_BLOCK),
            "row_block_rows": int(ROW_BLOCK),
            "dense_observation_gene_matrix": False,
        },
        "environment": _environment(),
    }
    return report


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Ticket #8 — five-fold simple baselines on the frozen 06.3 input.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT / STEP_REL / "06.3_input" / "output",
        help="ticket #3 frozen input directory",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=ROOT / STEP_REL,
        help="step container root; output goes to <out-root>/06.8_baselines/output",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="exact output directory (default: <out-root>/06.8_baselines/output)",
    )
    parser.add_argument(
        "--folds", type=int, nargs="+", default=None,
        help="folds to run (default: every frozen fold)",
    )
    parser.add_argument(
        "--baselines", nargs="+", choices=list(BASELINES), default=list(BASELINES),
        help="baselines to fit and score",
    )
    parser.add_argument(
        "--alpha", type=float, default=RIDGE_ALPHA,
        help="fixed ridge L2 strength (never tuned)",
    )
    parser.add_argument("--hash-cache", type=Path, default=None)
    parser.add_argument(
        "--table-cache", type=Path, default=None,
        help="derived CSV table cache (the pinned env has no parquet engine)",
    )
    parser.add_argument("--force", action="store_true", help="re-run existing fold/baseline outputs")
    parser.add_argument("--check", action="store_true", help="verify the frozen input, write nothing")
    parser.add_argument("--all", action="store_true", help="run every fold and write the tracked record")
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ns = build_parser().parse_args(argv)
    out_dir = Path(ns.out) if ns.out is not None else Path(ns.out_root) / "06.8_baselines" / "output"
    command = ["baselines.py", *argv]
    try:
        if ns.check:
            report = check_input(
                data_dir=ns.data_dir,
                folds=ns.folds,
                hash_cache=ns.hash_cache,
                table_cache=ns.table_cache,
            )
            print(json.dumps(report, indent=2, default=str))
            return 0
        if not ns.all:
            build_parser().error("one of --check / --all is required")
        run_all(
            data_dir=ns.data_dir,
            out_dir=out_dir,
            folds=ns.folds,
            baselines=ns.baselines,
            alpha=ns.alpha,
            hash_cache=ns.hash_cache,
            table_cache=ns.table_cache,
            force=ns.force,
            command=command,
        )
        return 0
    except (BaselineError, data_mod.DataError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
