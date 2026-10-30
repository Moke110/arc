"""Metric primitives for the Stage 2 `03_bench_models` contract (ticket 07).

All functions operate on float64 arrays unless stated otherwise.  The module is
import-safe: it has no side effects and no GPU dependency, so the build's
sanity cases and the offline tests can use it directly.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import average_precision_score


def true_set_indices(delta: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(up_idx, down_idx)`` for one group's true Δ.

    Top-*k* up and bottom-*k* down by the group-mean Δ.  Ties at the boundary
    are broken by corpus gene order (stable argsort), and a direction with
    fewer than *k* non-zero genes falls back to all non-zero genes.  The
    returned arrays are exactly *k* long on the current corpus.
    """

    delta = np.asarray(delta, dtype=np.float64)
    if delta.ndim != 1:
        raise ValueError(f"delta must be 1-D, got {delta.shape}")
    n = delta.size
    if not 1 <= k <= n:
        raise ValueError(f"k={k} outside 1..{n}")

    order_up = np.argsort(-delta, kind="stable")
    positive = delta > 0
    if int(positive.sum()) >= k:
        up = order_up[:k]
        if not bool((delta[up] > 0).all()):
            raise AssertionError("top-k up contains a non-positive gene")
    else:
        up = np.flatnonzero(positive)

    order_down = np.argsort(delta, kind="stable")
    negative = delta < 0
    if int(negative.sum()) >= k:
        down = order_down[:k]
        if not bool((delta[down] < 0).all()):
            raise AssertionError("bottom-k down contains a non-negative gene")
    else:
        down = np.flatnonzero(negative)

    if np.intersect1d(up, down).size:
        raise AssertionError("up and down true sets overlap")
    return up.astype(np.int64), down.astype(np.int64)


def direction_labels(n_genes: int, up: np.ndarray, down: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Build the two binary label vectors for one group."""

    y_up = np.zeros(n_genes, dtype=np.int8)
    y_down = np.zeros(n_genes, dtype=np.int8)
    y_up[up] = 1
    y_down[down] = 1
    return y_up, y_down


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(pearsonr(a, b).statistic)


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(spearmanr(a, b).statistic)


def r2_ss(true: np.ndarray, pred: np.ndarray) -> float:
    """``1 - SS_res / SS_tot`` about the true mean (never Pearson²)."""

    true = np.asarray(true, dtype=np.float64)
    pred = np.asarray(pred, dtype=np.float64)
    ss_res = float(np.sum((pred - true) ** 2))
    ss_tot = float(np.sum((true - true.mean()) ** 2))
    if ss_tot == 0.0:
        return math.nan
    return 1.0 - ss_res / ss_tot


def rmse(true: np.ndarray, pred: np.ndarray) -> float:
    true = np.asarray(true, dtype=np.float64)
    pred = np.asarray(pred, dtype=np.float64)
    return float(np.sqrt(np.mean((pred - true) ** 2)))


def nrmse(true: np.ndarray, pred: np.ndarray) -> float:
    sd = float(np.std(np.asarray(true, dtype=np.float64), ddof=0))
    if sd == 0.0:
        return math.nan
    return rmse(true, pred) / sd


def group_metrics(true: np.ndarray, pred: np.ndarray, k: int) -> dict[str, float]:
    """All per-group metrics at one *k*.

    ``auprc`` is the direction macro mean; ``auprc_up`` / ``auprc_down`` are
    persisted separately for audit.  A NaN direction is omitted from the mean;
    both-NaN yields NaN and the caller counts it.
    """

    true = np.asarray(true, dtype=np.float64)
    pred = np.asarray(pred, dtype=np.float64)
    up, down = true_set_indices(true, k)
    y_up, y_down = direction_labels(true.size, up, down)

    ap_up = float(average_precision_score(y_up, pred))
    ap_down = float(average_precision_score(y_down, -pred))
    vals = [v for v in (ap_up, ap_down) if not math.isnan(v)]
    ap = float(np.mean(vals)) if vals else math.nan

    return {
        "auprc": ap,
        "auprc_up": ap_up,
        "auprc_down": ap_down,
        "pcc": _pearson(true, pred),
        "spearman": _spearman(true, pred),
        "r2": r2_ss(true, pred),
        "rmse": rmse(true, pred),
        "nrmse": nrmse(true, pred),
    }
