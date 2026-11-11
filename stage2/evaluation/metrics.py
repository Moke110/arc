"""Canonical full-CM significance retrieval for every active Stage 2 step."""
from dataclasses import dataclass
import numpy as np
import pandas as pd

METRICS = ('AUPRC', 'AUPRC_positive', 'AUPRC_negative')
PR_FIELDS = ('recall', 'precision', 'contributing_groups', 'eligible_groups')


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
    """Signed covered truth and a full-axis positive count; zeros abstain."""
    delta = np.asarray(delta, dtype=np.float64)
    truth = np.asarray(significant, dtype=np.float64)
    if (delta.ndim != 1 or truth.shape != delta.shape or not np.isfinite(delta).all()
            or not np.isin(truth, (-1, 0, 1)).all() or positives < np.count_nonzero(truth)):
        raise ValueError('expected finite predictions, signed truth and full-axis positive count')
    if positives == 0:
        return None
    active = np.flatnonzero(delta != 0)
    ranked = active[np.argsort(-np.abs(delta[active]), kind='stable')]
    scores = np.abs(delta[ranked])
    if not len(ranked):
        empty = np.array([], dtype=float)
        return PR(0., np.array([0.]), np.array([0.]), empty, empty, empty, int(positives))
    end = np.r_[np.flatnonzero(np.diff(scores)), len(scores) - 1]
    correct = (truth[ranked] != 0) & (np.sign(delta[ranked]) == truth[ranked])
    tp = np.cumsum(correct, dtype=np.int64)[end]
    selected = end + 1
    recall, precision = tp / positives, tp / selected
    ap = float(np.dot(np.diff(np.r_[0., recall]), precision))
    return PR(ap, recall, precision, scores[end], selected, tp, int(positives))


def retrieval(prediction, signed_truth, covered_columns):
    """Score covered CM predictions against the complete measured CM truth."""
    prediction = np.asarray(prediction, dtype=float)
    truth = np.asarray(signed_truth)
    columns = np.asarray(covered_columns, dtype=int)
    if (len(columns) != len(prediction) or len(np.unique(columns)) != len(columns)
            or np.any(columns < 0) or np.any(columns >= len(truth))
            or not np.isin(truth, (-1, 0, 1)).all()):
        raise ValueError('invalid full-axis truth or prediction coverage')
    covered = truth[columns]
    result = {'AUPRC': average_precision(prediction, covered, np.count_nonzero(truth))}
    for name, sign in (('positive', 1), ('negative', -1)):
        keep = sign * prediction > 0
        labels = np.where(covered[keep] == sign, sign, 0)
        result[f'AUPRC_{name}'] = average_precision(prediction[keep], labels, np.count_nonzero(truth == sign))
    return result


class PooledPR:
    """Stream fixed-population threshold counts with constant memory per run."""
    def __init__(self):
        self.thresholds = np.r_[np.geomspace(1e4, 1e-12, 2048), 0.]
        self.hits = np.zeros(len(self.thresholds), dtype=np.int64)
        self.selected = np.zeros_like(self.hits)
        self.positives = self.groups = 0
        self.maximum = 0.
        self.maximum_hits = self.maximum_selected = 0

    def add(self, curve):
        if curve is None:
            return
        self.groups += 1
        self.positives += curve.positives
        if not len(curve.thresholds):
            return
        take = np.searchsorted(-curve.thresholds, -self.thresholds, side='right') - 1
        valid = take >= 0
        self.hits[valid] += curve.hits[take[valid]]
        self.selected[valid] += curve.selected[take[valid]]
        maximum = curve.thresholds[0]
        if maximum > self.maximum:
            self.maximum = maximum
            self.maximum_hits = self.maximum_selected = 0
        if maximum == self.maximum:
            self.maximum_hits += int(curve.hits[0])
            self.maximum_selected += int(curve.selected[0])

    def frame(self):
        if not self.groups:
            return pd.DataFrame(columns=PR_FIELDS)
        if not self.maximum_selected:
            return pd.DataFrame([[0., 0., self.groups, self.groups]], columns=PR_FIELDS)
        selected = np.r_[self.maximum_selected, self.selected]
        hits = np.r_[self.maximum_hits, self.hits]
        # Retain the actual highest threshold and every distinct sampled set.
        valid = selected > 0
        selected, hits = selected[valid], hits[valid]
        keep = np.r_[True, np.diff(selected) != 0]
        return pd.DataFrame({'recall': hits[keep] / self.positives,
                             'precision': hits[keep] / selected[keep],
                             'contributing_groups': self.groups, 'eligible_groups': self.groups}, columns=PR_FIELDS)


def pooled_pr(curves):
    result = PooledPR()
    for curve in curves:
        result.add(curve)
    return result.frame()
