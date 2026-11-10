"""Score one Round 1 group on its native axis and the common panel."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import rankdata


CONTINUOUS = ("pcc", "spearman", "rmse", "r2")
PR_FIELDS = ("recall", "precision", "contributing_groups", "eligible_groups")


@dataclass(frozen=True)
class PR:
    ap: float
    recall: np.ndarray
    precision: np.ndarray


def average_precision(delta, significant, positives):
    """Advance one threshold for all tied scores and retain only attained recall."""
    if positives == 0:
        return None
    scores = np.abs(delta)
    if len(scores) == 0:
        return PR(0., np.array([0.]), np.array([0.]))
    ranked = np.argsort(-scores, kind="stable")
    end = np.r_[np.flatnonzero(np.diff(scores[ranked])), len(ranked) - 1]
    tp = np.cumsum(significant[ranked], dtype=np.int64)[end]
    recall = tp / positives
    precision = tp / (end + 1)
    return PR(float(np.dot(np.diff(np.r_[0., recall]), precision)), recall, precision)


def macro_pr(curves):
    """Average precision at each recall only over groups reaching that recall."""
    if not curves:
        return pd.DataFrame(columns=PR_FIELDS)
    terminal = np.array([curve.recall[-1] for curve in curves])
    grid = np.unique(np.r_[np.linspace(0, 1, 101), terminal])
    grid = grid[grid <= terminal.max()]
    count = np.zeros(len(grid), dtype=np.int64)
    total = np.zeros(len(grid))
    for curve, limit in zip(curves, terminal, strict=True):
        supported = grid <= limit
        take = curve.recall.searchsorted(grid[supported], side="left")
        total[supported] += curve.precision[take]
        count[supported] += 1
    return pd.DataFrame({"recall": grid, "precision": total / count, "contributing_groups": count,
                         "eligible_groups": len(curves)}, columns=PR_FIELDS)


def correlation(pred, obs):
    if len(obs) < 2 or np.ptp(pred) == 0 or np.ptp(obs) == 0:
        return None
    a, b = pred - pred.mean(), obs - obs.mean()
    return float(np.clip((a @ b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))


def continuous(pred, obs):
    if not len(obs):
        return dict.fromkeys(CONTINUOUS)
    error = np.square(pred - obs).sum(dtype=np.float64)
    variance = np.square(obs - obs.mean()).sum(dtype=np.float64)
    return {"pcc": correlation(pred, obs), "spearman": correlation(rankdata(pred), rankdata(obs)),
            "rmse": float(np.sqrt(error / len(obs))), "r2": float(1 - error / variance) if variance else None}


def project_delta(delta, native_control, columns, shared_control):
    log_expression = np.maximum(delta[:, columns].astype(np.float64) + native_control[:, columns], 0.)
    if not log_expression.shape[1]:
        return log_expression
    # expm1(log_expression), divided by a common scale to avoid overflow.
    linear = np.exp(log_expression - log_expression.max(axis=1, keepdims=True)) * -np.expm1(-log_expression)
    sums = linear.sum(axis=1, keepdims=True)
    np.divide(linear, sums, out=linear, where=sums > 0)
    return np.log1p(10_000 * linear) - shared_control


def paired_promotions(rows):
    scores = {row["model"]: row["auprc_full"] for row in rows}
    decisions = []
    for model, value in scores.items():
        if model.startswith("XPert_G"):
            control = "XPert_base"
        elif model.startswith("XPert_MLP"):
            control = "XPert_G18164"
        else:
            continue
        control_value = scores.get(control)
        if value is None or control_value is None:
            decision = "unavailable"
        else:
            decision = "promoted" if value > control_value else "not_promoted"
        decisions.append({"model": model, "control": control, "decision": decision,
                          "model_auprc_full": value, "control_auprc_full": control_value})
    return decisions
