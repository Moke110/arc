"""Direction AUPRC@k and the auxiliary Δ regression metrics.

**Same-source copy — do not diverge** for the metric *definitions*: the
per-direction ranking-label construction (``true_set_indices`` /
``direction_labels`` / ``group_metrics``) is copied from
``stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/metrics.py``,
which implements the frozen benchmark contract:

* the true Δ vector of a condition is ranked; the top ``k`` genes are the "up"
  labels and the bottom ``k`` the "down" labels, ties broken stably by gene
  order (a direction with fewer than ``k`` non-zero genes falls back to all
  non-zero genes in that direction);
* ``AP_up = average_precision_score(labels_up, pred)`` and
  ``AP_down = average_precision_score(labels_down, -pred)``; ``auprc`` is the
  mean over the non-NaN directions;
* ``r2`` is ``1 - SS_res/SS_tot`` about the true mean, never Pearson².

What this module adds is a **vectorised, tie-aware** implementation
(``average_precision_matrix``) so the per-epoch validation score of both stages
costs a chunked argsort + cumulative sum instead of a Python loop over
conditions.  Its AP values agree with ``sklearn.metrics.average_precision_score``
including exact ties (checked by ``tests/test_metrics.py``).

Aggregation order follows the spec: within a condition, replicate predictions
are averaged **first** (``aggregate_rows`` / the caller), then the per-condition
direction APs are computed, then macro-averaged over conditions.
"""

from __future__ import annotations

import math
import warnings
from typing import Any, Iterable, Mapping

import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import average_precision_score

__all__ = [
    "K_LIST",
    "true_set_indices",
    "direction_labels",
    "group_metrics",
    "true_set_masks",
    "average_precision_matrix",
    "direction_auprc",
    "direction_auprc_matrix",
    "direction_auprc_conditions",
    "regression_metrics_matrix",
    "aggregate_rows",
    "aggregate_metrics",
    "macro",
    "DEFAULT_CHUNK",
]

#: k values persisted for audit; k=200 is the headline (``config.metric.k``).
K_LIST = (100, 200, 500)
DEFAULT_CHUNK = 256


# --------------------------------------------------------------------------- #
# Frozen per-condition label construction (copied from 03.1 metrics.py)
# --------------------------------------------------------------------------- #
def true_set_indices(delta: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(up_idx, down_idx)`` for one condition's true Δ."""
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


def direction_labels(
    n_genes: int, up: np.ndarray, down: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Build the two binary label vectors for one condition."""
    y_up = np.zeros(n_genes, dtype=np.int8)
    y_down = np.zeros(n_genes, dtype=np.int8)
    y_up[up] = 1
    y_down[down] = 1
    return y_up, y_down


def true_set_masks(delta: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised :func:`true_set_indices` over a ``(conditions, genes)`` matrix.

    Returns ``(mask_up, mask_down)`` int8 matrices.  The fallback (fewer than
    ``k`` non-zero genes in a direction) is handled row by row; it is rare and
    kept out of the vectorised hot path.
    """
    delta = np.asarray(delta, dtype=np.float64)
    if delta.ndim != 2:
        raise ValueError(f"delta must be 2-D (conditions, genes), got {delta.shape}")
    n_cond, n_genes = delta.shape
    if not 1 <= k <= n_genes:
        raise ValueError(f"k={k} outside 1..{n_genes}")

    n_pos = (delta > 0).sum(axis=1)
    n_neg = (delta < 0).sum(axis=1)
    up = np.zeros((n_cond, n_genes), dtype=np.int8)
    down = np.zeros((n_cond, n_genes), dtype=np.int8)

    order_up = np.argsort(-delta, axis=1, kind="stable")
    full_up = np.flatnonzero(n_pos >= k)
    if full_up.size:
        up[full_up[:, None], order_up[full_up, :k]] = 1
    for row in np.flatnonzero(n_pos < k):
        up[row, delta[row] > 0] = 1

    order_down = np.argsort(delta, axis=1, kind="stable")
    full_down = np.flatnonzero(n_neg >= k)
    if full_down.size:
        down[full_down[:, None], order_down[full_down, :k]] = 1
    for row in np.flatnonzero(n_neg < k):
        down[row, delta[row] < 0] = 1

    if np.any(up & down):
        raise AssertionError("up and down true sets overlap")
    return up, down


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
    return float(np.sqrt(np.mean((np.asarray(pred, np.float64) - np.asarray(true, np.float64)) ** 2)))


def nrmse(true: np.ndarray, pred: np.ndarray) -> float:
    sd = float(np.std(np.asarray(true, dtype=np.float64), ddof=0))
    if sd == 0.0:
        return math.nan
    return rmse(true, pred) / sd


def group_metrics(true: np.ndarray, pred: np.ndarray, k: int) -> dict[str, float]:
    """All per-condition metrics at one *k* (single condition; audit path)."""
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


# --------------------------------------------------------------------------- #
# Vectorised AP (tie-aware, identical to sklearn average_precision_score)
# --------------------------------------------------------------------------- #
def average_precision_matrix(scores: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Per-row average precision for a ``(conditions, genes)`` score matrix.

    Implements the same definition as
    ``sklearn.metrics.average_precision_score``: rank genes by descending score,
    collapse exact ties into one threshold, and sum
    ``sum_i (R_i - R_{i-1}) * P_i`` where ``P_i`` is the precision at the end of
    tie group *i*.  Rows with no positive label return ``NaN``.
    """
    scores = np.asarray(scores, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int8)
    if scores.shape != labels.shape or scores.ndim != 2:
        raise ValueError("scores and labels must be 2-D and the same shape")
    n_cond, n_genes = scores.shape
    n_pos = labels.sum(axis=1).astype(np.float64)

    order = np.argsort(-scores, axis=1, kind="stable")
    sorted_scores = np.take_along_axis(scores, order, axis=1)
    sorted_labels = np.take_along_axis(labels.astype(np.int16), order, axis=1)
    cum_pos = np.cumsum(sorted_labels, axis=1, dtype=np.int64)

    is_end = np.ones((n_cond, n_genes), dtype=bool)
    if n_genes > 1:
        is_end[:, :-1] = sorted_scores[:, :-1] != sorted_scores[:, 1:]
    ranks = np.arange(1, n_genes + 1, dtype=np.float64)
    precision = cum_pos / ranks

    # Forward-fill "the end index of my tie group": ends are increasing
    # positions, so the group end at or after j is the minimum of them.
    end_index = np.where(is_end, np.arange(n_genes), n_genes)
    end_index = np.minimum.accumulate(end_index[:, ::-1], axis=1)[:, ::-1]
    precision_at_end = np.take_along_axis(precision, end_index, axis=1)

    contributions = precision_at_end * sorted_labels.astype(np.float64)
    ap = contributions.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        ap = np.where(n_pos > 0, ap / n_pos, np.nan)
    return ap


# --------------------------------------------------------------------------- #
# Direction AUPRC over a condition matrix
# --------------------------------------------------------------------------- #
def direction_auprc_matrix(
    true_delta: np.ndarray,
    pred_delta: np.ndarray,
    k: int,
    chunk: int = DEFAULT_CHUNK,
) -> dict[str, np.ndarray]:
    """Per-condition direction AUPRC for ``(conditions, genes)`` Δ matrices.

    Returns arrays of length ``conditions``: ``auprc`` (mean over non-NaN
    directions), ``auprc_up``, ``auprc_down``, plus the true-set sizes.
    """
    true_delta = np.asarray(true_delta, dtype=np.float64)
    pred_delta = np.asarray(pred_delta, dtype=np.float64)
    if true_delta.shape != pred_delta.shape or true_delta.ndim != 2:
        raise ValueError("true_delta and pred_delta must be 2-D and the same shape")
    n_cond = true_delta.shape[0]
    auprc = np.full(n_cond, np.nan, dtype=np.float64)
    auprc_up = np.full(n_cond, np.nan, dtype=np.float64)
    auprc_down = np.full(n_cond, np.nan, dtype=np.float64)
    n_up = np.zeros(n_cond, dtype=np.int64)
    n_down = np.zeros(n_cond, dtype=np.int64)

    for start in range(0, n_cond, chunk):
        stop = min(start + chunk, n_cond)
        t = true_delta[start:stop]
        p = pred_delta[start:stop]
        mask_up, mask_down = true_set_masks(t, k)
        ap_up = average_precision_matrix(p, mask_up)
        ap_down = average_precision_matrix(-p, mask_down)
        auprc_up[start:stop] = ap_up
        auprc_down[start:stop] = ap_down
        stacked = np.vstack([ap_up, ap_down])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            auprc[start:stop] = np.nanmean(stacked, axis=0)
        n_up[start:stop] = mask_up.sum(axis=1)
        n_down[start:stop] = mask_down.sum(axis=1)

    return {
        "auprc": auprc,
        "auprc_up": auprc_up,
        "auprc_down": auprc_down,
        "n_up": n_up,
        "n_down": n_down,
    }


def direction_auprc(
    true_delta: np.ndarray, pred_delta: np.ndarray, k: int
) -> dict[str, np.ndarray]:
    """Alias of :func:`direction_auprc_matrix` (explicit single-vector entry point)."""
    true_delta = np.asarray(true_delta, dtype=np.float64)
    pred_delta = np.asarray(pred_delta, dtype=np.float64)
    squeeze = true_delta.ndim == 1
    if squeeze:
        true_delta = true_delta[None, :]
        pred_delta = pred_delta[None, :]
    out = direction_auprc_matrix(true_delta, pred_delta, k)
    if squeeze:
        return {key: value[0] for key, value in out.items()}
    return out


def direction_auprc_conditions(
    conditions: Iterable[str],
    true_delta: np.ndarray,
    pred_delta: np.ndarray,
    k: int,
    chunk: int = DEFAULT_CHUNK,
) -> list[dict[str, Any]]:
    """Per-condition records (the persisted ``per_condition.csv`` rows)."""
    names = [str(c) for c in conditions]
    if len(names) != np.asarray(true_delta).shape[0]:
        raise ValueError("condition count does not match the Δ matrix rows")
    metrics = direction_auprc_matrix(true_delta, pred_delta, k, chunk=chunk)
    records = []
    for i, name in enumerate(names):
        records.append(
            {
                "condition": name,
                "auprc": float(metrics["auprc"][i]),
                "auprc_up": float(metrics["auprc_up"][i]),
                "auprc_down": float(metrics["auprc_down"][i]),
                "n_up": int(metrics["n_up"][i]),
                "n_down": int(metrics["n_down"][i]),
            }
        )
    return records


# --------------------------------------------------------------------------- #
# Auxiliary Δ regression metrics (vectorised over conditions)
# --------------------------------------------------------------------------- #
def _row_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean(axis=1, keepdims=True)
    denom = np.sqrt((a * a).sum(axis=1) * (b * b).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (a * b).sum(axis=1) / denom
    out[denom == 0.0] = np.nan
    return out


def _row_rank(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, axis=1, kind="stable")
    ranks = np.empty(values.shape, dtype=np.float64)
    np.put_along_axis(
        ranks, order, np.broadcast_to(np.arange(values.shape[1], dtype=np.float64), values.shape), axis=1
    )
    return ranks


def regression_metrics_matrix(
    true_delta: np.ndarray, pred_delta: np.ndarray, chunk: int = DEFAULT_CHUNK
) -> dict[str, np.ndarray]:
    """Δ RMSE / Pearson / Spearman per condition (vectorised, chunked)."""
    true_delta = np.asarray(true_delta, dtype=np.float64)
    pred_delta = np.asarray(pred_delta, dtype=np.float64)
    if true_delta.shape != pred_delta.shape or true_delta.ndim != 2:
        raise ValueError("true_delta and pred_delta must be 2-D and the same shape")
    n_cond = true_delta.shape[0]
    out = {
        "rmse": np.full(n_cond, np.nan),
        "pcc": np.full(n_cond, np.nan),
        "spearman": np.full(n_cond, np.nan),
        "pred_constant": np.zeros(n_cond, dtype=bool),
        "true_constant": np.zeros(n_cond, dtype=bool),
    }
    for start in range(0, n_cond, chunk):
        stop = min(start + chunk, n_cond)
        t = true_delta[start:stop]
        p = pred_delta[start:stop]
        diff = p - t
        out["rmse"][start:stop] = np.sqrt((diff * diff).mean(axis=1))
        pred_constant = np.all(p == p[:, :1], axis=1)
        true_constant = np.all(t == t[:, :1], axis=1)
        out["pred_constant"][start:stop] = pred_constant
        out["true_constant"][start:stop] = true_constant
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            pcc = _row_pearson(t, p)
            pcc[pred_constant | true_constant] = np.nan
            sp = _row_pearson(_row_rank(t), _row_rank(p))
            sp[pred_constant | true_constant] = np.nan
        out["pcc"][start:stop] = pcc
        out["spearman"][start:stop] = sp
    return out


# --------------------------------------------------------------------------- #
# Aggregation
# --------------------------------------------------------------------------- #
def aggregate_rows(values: np.ndarray, group_index: np.ndarray, n_groups: int) -> np.ndarray:
    """Mean of rows within each group — the "aggregate replicates first" step.

    ``values`` is ``(rows, ...)`` and ``group_index`` maps each row to its
    condition.  This is the numpy twin of 03.1's ``group_mean_predictions`` and
    is applied to Δ (predictions and truth) *before* any ranking metric.
    """
    values = np.asarray(values)
    group_index = np.asarray(group_index)
    if values.shape[0] != group_index.size:
        raise ValueError("row count does not match the group index")
    if group_index.size and (group_index.min() < 0 or group_index.max() >= n_groups):
        raise ValueError("group index out of range")
    flat = values.reshape(values.shape[0], -1)
    sums = np.zeros((n_groups, flat.shape[1]), dtype=np.float64)
    counts = np.zeros(n_groups, dtype=np.int64)
    np.add.at(sums, group_index, flat.astype(np.float64))
    np.add.at(counts, group_index, 1)
    if np.any(counts == 0):
        raise ValueError("empty condition group")
    return (sums / counts[:, None]).reshape((n_groups, *values.shape[1:]))


def macro(values: Any) -> dict[str, float]:
    """Macro aggregate of a metric vector: mean/sd/min/nan, NaN omitted."""
    arr = np.asarray(values, dtype=np.float64).reshape(-1)
    finite = arr[np.isfinite(arr)]
    n_nan = int(arr.size - finite.size)
    if finite.size == 0:
        return {"mean": math.nan, "sd": math.nan, "min": math.nan, "nan": n_nan}
    sd = float(np.std(finite, ddof=1)) if finite.size > 1 else math.nan
    return {
        "mean": float(np.mean(finite)),
        "sd": sd,
        "min": float(np.min(finite)),
        "nan": n_nan,
    }


def aggregate_metrics(metrics: Mapping[str, np.ndarray]) -> dict[str, dict[str, float]]:
    """Macro-average every metric array returned by the scorers."""
    out: dict[str, dict[str, float]] = {}
    for key, values in metrics.items():
        if key in {"n_up", "n_down", "pred_constant", "true_constant"}:
            continue
        out[key] = macro(values)
    return out
