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
    thresholds: np.ndarray
    selected: np.ndarray
    hits: np.ndarray
    positives: int


def average_precision(delta, significant, positives):
    """Signed truth: -1/down, +1/up, 0/nonsignificant; zero predictions abstain."""
    if positives == 0:
        return None
    active = np.flatnonzero(np.abs(delta) > 0)
    ranked = active[np.argsort(-np.abs(delta[active]), kind="stable")]
    scores = np.abs(delta[ranked])
    if not len(ranked):
        empty = np.array([], dtype=float)
        return PR(0., np.array([0.]), np.array([0.]), empty, empty, empty, int(positives))
    end = np.r_[np.flatnonzero(np.diff(scores)), len(ranked) - 1]
    truth = np.asarray(significant, dtype=np.float64)[ranked]
    correct = (truth != 0) & (np.sign(delta[ranked]) == np.sign(truth))
    tp = np.cumsum(correct, dtype=np.int64)[end]
    selected = end + 1
    recall = tp / positives
    precision = tp / selected
    return PR(float(np.dot(np.diff(np.r_[0., recall]), precision)), recall, precision,
              scores[end], selected, tp, int(positives))


def pooled_pr(curves):
    """Pool threshold events across a fixed group population (micro PR).

    Macro AP remains the mean of exact per-group AP; it is not the area
    of this pooled display curve.
    """
    if not curves:
        return pd.DataFrame(columns=PR_FIELDS)
    count = len(curves)
    active = [c for c in curves if len(c.thresholds)]
    if not active:
        return pd.DataFrame([[0., 0., count, count]], columns=PR_FIELDS)
    thresholds = np.concatenate([c.thresholds for c in active])
    selected = np.concatenate([np.diff(np.r_[0, c.selected]) for c in active])
    hits = np.concatenate([np.diff(np.r_[0, c.hits]) for c in active])
    order = np.argsort(-thresholds, kind="stable")
    end = np.r_[np.flatnonzero(np.diff(thresholds[order])), len(order) - 1]
    tp = np.cumsum(hits[order], dtype=np.int64)[end]
    retrieved = np.cumsum(selected[order], dtype=np.int64)[end]
    recall = tp / sum(c.positives for c in curves)
    precision = tp / retrieved
    # Bound published curve size, retaining actual thresholds and the endpoint.
    keep = np.unique(np.r_[0, np.linspace(0, len(end) - 1, min(2001, len(end))).astype(int), len(end) - 1])
    return pd.DataFrame({"recall": recall[keep], "precision": precision[keep],
                         "contributing_groups": count, "eligible_groups": count}, columns=PR_FIELDS)


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
