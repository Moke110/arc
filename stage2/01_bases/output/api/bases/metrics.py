"""Indexed retrieval and continuous metrics owned by the standalone bases."""
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
    original_columns = np.asarray(covered_columns)
    if (prediction.ndim != 1 or truth.ndim != 1 or original_columns.ndim != 1
            or (original_columns.size and original_columns.dtype.kind not in 'iu')):
        raise ValueError('prediction coverage requires integer gene indices on a one-dimensional axis')
    columns = original_columns.astype(np.int64)
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


def continuous(truth, prediction):
    """PCC, Spearman, RMSE and R2 on one declared covered numeric axis."""
    from scipy.stats import rankdata
    truth, prediction = np.asarray(truth, dtype=float), np.asarray(prediction, dtype=float)
    if (truth.ndim != 1 or truth.shape != prediction.shape or not len(truth)
            or not np.isfinite(truth).all() or not np.isfinite(prediction).all()):
        raise ValueError('continuous truth/prediction must share a finite nonempty axis')
    a, b = truth - truth.mean(), prediction - prediction.mean()
    variance, other = float(a @ a), float(b @ b)
    result = {'PCC': None, 'Spearman': None}
    if variance and other:
        result['PCC'] = float((a @ b) / np.sqrt(variance * other))
        result['Spearman'] = float(np.corrcoef(rankdata(truth), rankdata(prediction))[0, 1])
    error = float((truth - prediction) @ (truth - prediction))
    result['RMSE'] = float(np.sqrt(error / len(truth)))
    result['R2'] = 1 - error / variance if variance else None
    return result


def summarize_runs(summaries):
    """Equal-run mean of eligible group-macro AP values, never pooled rows."""
    summaries = list(summaries)
    result = {}
    for name in METRICS:
        values = [x[name] for x in summaries if x.get(name) is not None]
        if not np.isfinite(values).all():
            raise ValueError('nonfinite run metric')
        result[name] = float(np.mean(values)) if values else None
    return result


def aggregate_predictions(predictions, samples):
    """Yield means of well deltas for each full parent-dose group."""
    required = {'sample_id', 'drug_parent', 'dose_um'}
    if not required <= set(samples) or samples.sample_id.isna().any() or samples.sample_id.duplicated().any():
        raise ValueError('unique sample identities and parent-dose metadata are required')
    if (predictions.ndim != 2 or predictions.shape[0] != len(samples)
            or samples.drug_parent.isna().any() or samples.drug_parent.astype(str).eq('').any()
            or not np.isfinite(samples.dose_um).all() or samples.dose_um.le(0).any()):
        raise ValueError('prediction rows or normalized parent-dose metadata differ')
    rows = samples.reset_index(drop=True)
    dose_keys = rows.dose_um.map(lambda dose:format(float(dose),'.12g'))
    for (parent, dose), part in rows.groupby([rows.drug_parent,dose_keys], sort=True):
        values = np.asarray(predictions[part.index.to_numpy()], dtype=np.float64)
        if not np.isfinite(values).all():
            raise ValueError('nonfinite well prediction')
        yield str(parent), float(dose), part.sample_id.astype(str).tolist(), values.mean(axis=0)


def score_groups(predictions, samples, truth, covered_columns, *, source='tenv'):
    """Partition-authorized measured truth; predicted group means are not filtered.

    truth implements PartitionTruth.truth(parent, dose, exact_sample_ids).
    Untestable truth returns None; tested-empty arrays yield undefined PR.
    """
    if source != 'tenv':
        raise ValueError('measured significance retrieval is tenv-only')
    pools = {name: PooledPR() for name in METRICS}
    rows = []
    for parent, dose, members, prediction in aggregate_predictions(predictions, samples):
        signed = truth.truth(parent, dose, members)
        curves = {name: None for name in METRICS} if signed is None else retrieval(prediction, signed, covered_columns)
        row = {'drug_parent': parent, 'dose_um': dose, 'n_wells': len(members),
               'status': 'untestable' if signed is None else 'tested'}
        for name, curve in curves.items():
            row[name] = None if curve is None else curve.ap
            pools[name].add(curve)
        rows.append(row)
    frame = pd.DataFrame(rows, columns=['drug_parent', 'dose_um', 'n_wells', 'status', *METRICS])
    summary = {'groups': len(rows), 'eligible_groups': {}}
    for name in METRICS:
        values = frame[name].dropna().to_numpy(dtype=float)
        summary[name] = float(values.mean()) if len(values) else None
        summary['eligible_groups'][name] = len(values)
    return {'groups': frame, 'summary': summary, 'curves': {name: pool.frame() for name, pool in pools.items()}}
