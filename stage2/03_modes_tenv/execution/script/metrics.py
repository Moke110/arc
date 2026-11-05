"""Tie-aware directional retrieval and continuous group metrics (ADR 0042)."""
import numpy as np
from scipy.stats import rankdata

METRICS = ('AUPRC_positive@100', 'AUPRC_negative@100', 'Recall_positive@P80',
           'Recall_negative@P80', 'PCC', 'Spearman', 'RMSE', 'R2')
RECALL_GRID = np.linspace(0, 1, 101)
CURVE_CONVENTION = ('Nonempty tied-score thresholds only; at grid recall r use precision '
                    'of first threshold attaining r (right step); no synthetic (0,1) endpoint. '
                    'Macro-average eligible groups within run before run mean/min/max.')


def directional_curve(labels, scores):
    """Return AP, attainable P80 recall and precision on RECALL_GRID.

    Tied scores enter together. AP is the sum of precision times recall increments
    over actual score thresholds, equivalent to non-interpolated average precision.
    """
    labels = np.asarray(labels, bool)
    scores = np.asarray(scores, dtype=np.float64)
    if labels.shape != scores.shape or labels.ndim != 1 or not np.isfinite(scores).all():
        raise ValueError('directional labels/scores require matching finite vectors')
    total = int(labels.sum())
    if not total:
        return np.nan, np.nan, np.full(len(RECALL_GRID), np.nan)
    order = np.argsort(-scores, kind='stable')
    ends = np.r_[np.flatnonzero(np.diff(scores[order]) != 0), len(scores) - 1]
    tp = np.cumsum(labels[order], dtype=np.int64)[ends]
    recall = tp / total
    precision = tp / (ends + 1)
    ap = float(np.dot(np.diff(np.r_[0., recall]), precision))
    r80 = float(recall[precision >= .8].max(initial=0))
    indices = np.searchsorted(recall, RECALL_GRID, side='left')
    curve = precision[np.minimum(indices, len(precision) - 1)]
    return ap, r80, curve


def truth_labels(truth):
    truth = np.asarray(truth)
    labels = []
    for sign in (1, -1):
        axis = sign * truth
        eligible = np.flatnonzero(axis > 0)
        selected = eligible[np.argsort(-axis[eligible], kind='stable')[:100]]
        label = np.zeros(len(truth), dtype=bool)
        label[selected] = True
        labels.append(label)
    return labels


def score_group(truth, prediction, labels=None):
    truth = np.asarray(truth, dtype=np.float64)
    prediction = np.asarray(prediction, dtype=np.float64)
    if truth.shape != prediction.shape or truth.ndim != 1 or not np.isfinite(truth).all() or not np.isfinite(prediction).all():
        raise ValueError('truth/prediction require matching finite vectors')
    labels = truth_labels(truth) if labels is None else labels
    positive = directional_curve(labels[0], prediction)
    negative = directional_curve(labels[1], -prediction)
    values = dict(zip(METRICS[:4], (positive[0], negative[0], positive[1], negative[1])))
    reasons = {}
    for direction, result in (('positive', positive), ('negative', negative)):
        if np.isnan(result[0]):
            reasons[f'AUPRC_{direction}@100'] = reasons[f'Recall_{direction}@P80'] = 'empty_direction'
    centered = truth - truth.mean()
    other = prediction - prediction.mean()
    variance = float(np.dot(centered, centered))
    predicted_variance = float(np.dot(other, other))
    if variance and predicted_variance:
        values['PCC'] = float(np.dot(centered, other) / np.sqrt(variance * predicted_variance))
        x, y = rankdata(truth), rankdata(prediction)
        values['Spearman'] = float(np.corrcoef(x, y)[0, 1])
    else:
        values['PCC'] = values['Spearman'] = np.nan
        reason = 'constant_truth' if not variance else 'constant_prediction'
        reasons['PCC'] = reasons['Spearman'] = reason
    squared_error = float(np.dot(truth - prediction, truth - prediction))
    values['RMSE'] = float(np.sqrt(squared_error / len(truth)))
    values['R2'] = 1. - squared_error / variance if variance else np.nan
    if not variance:
        reasons['R2'] = 'zero_variance_truth'
    return values, reasons, np.stack((positive[2], negative[2]))
