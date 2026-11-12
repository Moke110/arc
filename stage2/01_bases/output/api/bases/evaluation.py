"""Stage 2 evaluation semantics shared by every model-comparison step.

Owns the evaluator truth gate, the per-well direct metrics, the predicted
absolute rank table, the four group-level indirect metrics and run-level
aggregation. Every function is pure over arrays plus identity metadata: no file
I/O, no cached state. Undefined values are returned as NaN together with a
reason; nothing is zero-filled. The retired AUPRC/PR implementation in
``bases.metrics`` is deliberately not imported or modified here (ADR 0063).
"""
from dataclasses import dataclass
import numpy as np

GENE_COUNT = 18164
K_MAX = 100
DIRECT_METRICS = ('PCC', 'Spearman', 'RMSE', 'R2')
INDIRECT_METRICS = ('Precision@k', 'Spearman@k', 'max_k@P80S40', 'Recall@max_k')

PADJ_THRESHOLD = 0.05
EFFECT_THRESHOLD = 2.0
K_MIN_FOR_SPEARMAN = 3
PRECISION_TARGET = 0.8
SPEARMAN_TARGET = 0.4

# Undefined reasons; never used to fill a value with zero.
EMPTY_GATE = 'empty_gate'
K_BELOW_THREE = 'k_below_3'
NO_QUALIFYING_K = 'no_qualifying_k'
CONSTANT_PREDICTION = 'constant_prediction'
CONSTANT_TRUTH = 'constant_truth'
NO_DEFINED_UNITS = 'no_defined_units'
NO_DEFINED_GROUPS = 'no_defined_groups'


def _vector(value, label):
    array = np.asarray(value)
    if array.ndim != 1 or not len(array):
        raise ValueError(f'{label} must be a nonempty one-dimensional axis')
    return array


def _float_vector(value, label):
    array = _vector(value, label).astype(np.float64)
    if not np.isfinite(array).all():
        raise ValueError(f'{label} must be finite')
    return array


def _shared_axis(**arrays):
    values = {label: _vector(value, label) for label, value in arrays.items()}
    lengths = {len(value) for value in values.values()}
    if len(lengths) != 1:
        raise ValueError('evaluation arrays must share one gene axis')
    return values


@dataclass(frozen=True)
class GatedTruth:
    """Full-axis signed gate: +1/-1 for gated truth genes, 0 elsewhere."""
    signed: np.ndarray
    indices: np.ndarray

    @property
    def count(self):
        return int(len(self.indices))


def truth_gate(padj, log2_fold_change, delta_mean):
    """Gated truth: ``padj < 0.05`` and ``|log2FC| >= 2``, direction from Δ_mean."""
    arrays = _shared_axis(padj=padj, log2_fold_change=log2_fold_change, delta_mean=delta_mean)
    padj = arrays['padj'].astype(np.float64)
    effect = arrays['log2_fold_change'].astype(np.float64)
    delta = arrays['delta_mean'].astype(np.float64)
    finite = np.isfinite(padj)
    if np.any(np.isinf(padj)) or np.any((padj[finite] < 0) | (padj[finite] > 1)):
        raise ValueError('invalid adjusted p-values')
    selected = padj < PADJ_THRESHOLD
    if np.any(selected & (~np.isfinite(effect) | (effect == 0))):
        raise ValueError('significant gene has invalid direction')
    gated = selected & (np.abs(effect) >= EFFECT_THRESHOLD)
    if np.any(gated & (~np.isfinite(delta) | (delta == 0))):
        raise ValueError('gated truth gene has no measured delta direction')
    signed = np.zeros(len(padj), dtype=np.int8)
    indices = np.flatnonzero(gated).astype(np.int64)
    signed[indices] = np.sign(delta[indices]).astype(np.int8)
    return GatedTruth(signed=signed, indices=indices)


@dataclass(frozen=True)
class DirectMetrics:
    """One test well's direct metrics plus reasons for every undefined value."""
    values: dict
    reasons: dict


def direct_metrics(truth, prediction):
    """PCC, Spearman, RMSE and R2 on one well over the full gene axis.

    A constant truth or prediction leaves PCC, Spearman and R2 undefined: NaN
    with ``constant_truth``/``constant_prediction``. RMSE stays defined.
    """
    arrays = _shared_axis(truth=truth, prediction=prediction)
    truth = arrays['truth'].astype(np.float64)
    prediction = arrays['prediction'].astype(np.float64)
    if not np.isfinite(truth).all() or not np.isfinite(prediction).all():
        raise ValueError('direct metric truth and prediction must be finite')
    from scipy.stats import rankdata

    centered_truth, centered_prediction = truth - truth.mean(), prediction - prediction.mean()
    truth_variance, prediction_variance = float(centered_truth @ centered_truth), float(
        centered_prediction @ centered_prediction)
    values = {name: float('nan') for name in DIRECT_METRICS}
    reasons = {}
    if not truth_variance:
        reasons.update({name: CONSTANT_TRUTH for name in ('PCC', 'Spearman', 'R2')})
    elif not prediction_variance:
        reasons.update({name: CONSTANT_PREDICTION for name in ('PCC', 'Spearman', 'R2')})
    else:
        values['PCC'] = float(centered_truth @ centered_prediction) / np.sqrt(
            truth_variance * prediction_variance)
        values['Spearman'] = float(np.corrcoef(rankdata(truth), rankdata(prediction))[0, 1])
    error = float((truth - prediction) @ (truth - prediction))
    values['RMSE'] = float(np.sqrt(error / len(truth)))
    if truth_variance and prediction_variance:
        values['R2'] = 1 - error / truth_variance
    return DirectMetrics(values=values, reasons=reasons)


@dataclass(frozen=True)
class PredictedRankTable:
    """Full-axis absolute predicted rank table for one group."""
    gene_index: np.ndarray
    rank: np.ndarray
    direction: np.ndarray
    absolute_delta: np.ndarray

    def __len__(self):
        return int(len(self.gene_index))


def predicted_rank_table(delta_mean):
    """Order the full gene axis by ``|Δ̂_mean|`` descending.

    Ties in absolute predicted Δ are ordered by gene index so the table is
    deterministic, and every gene of the axis appears exactly once.
    """
    delta = _float_vector(delta_mean, 'delta_mean')
    gene_index = np.arange(len(delta), dtype=np.int64)
    absolute = np.abs(delta)
    order = np.lexsort((gene_index, -absolute))
    return PredictedRankTable(
        gene_index=gene_index[order],
        rank=np.arange(1, len(delta) + 1, dtype=np.int64),
        direction=np.sign(delta[order]).astype(np.int8),
        absolute_delta=absolute[order],
    )


@dataclass(frozen=True)
class GroupIndirectMetrics:
    """One group's indirect curves, qualifying cut-off and recall, with reasons.

    ``matched`` is the number of direction-matched gated truth genes on the whole
    axis, i.e. the Recall@max_k denominator.
    """
    k: np.ndarray
    precision: np.ndarray
    spearman: np.ndarray
    max_k: float
    recall: float
    matched: int
    reasons: dict


def indirect_metrics(table, gated, delta_mean, *, k_max=K_MAX):
    """Precision@k, Spearman@k, max_k@P80S40 and Recall@max_k for one group.

    ``table`` is the group's predicted absolute rank table, ``gated`` its truth
    gate and ``delta_mean`` the group's measured Δ_mean. ``k`` runs over
    1..``k_max``; every undefined value is NaN with a reason. Spearman@k is the
    rank correlation between signed predicted Δ̂_mean and measured Δ_mean over
    the top-k set (average ranks resolve ties).
    """
    from scipy.stats import rankdata

    delta = _float_vector(delta_mean, 'delta_mean')
    signed_truth = np.asarray(gated.signed)
    if (signed_truth.ndim != 1 or len(signed_truth) != len(delta) or len(table) != len(delta)):
        raise ValueError('indirect metrics require one shared gene axis')
    if not isinstance(k_max, (int, np.integer)) or k_max < 1:
        raise ValueError('k_max must be a positive integer')
    if k_max > len(delta):
        raise ValueError('k_max exceeds the gene axis')
    k_values = np.arange(1, k_max + 1, dtype=np.int64)
    precision = np.full(k_max, np.nan)
    spearman = np.full(k_max, np.nan)
    max_k = float('nan')
    recall = float('nan')
    reasons = {}
    ordered_truth = signed_truth[table.gene_index]
    ordered_delta = delta[table.gene_index]
    ordered_direction = np.asarray(table.direction)
    ordered_prediction = ordered_direction.astype(np.float64) * np.asarray(table.absolute_delta)
    matched = (ordered_truth != 0) & (ordered_direction == ordered_truth)
    matched_total = int(np.count_nonzero(matched))
    if not gated.count:
        reasons.update({name: EMPTY_GATE for name in INDIRECT_METRICS})
        return GroupIndirectMetrics(k_values, precision, spearman, max_k, recall, 0, reasons)
    cumulative = np.cumsum(matched, dtype=np.int64)[:k_max]
    precision = cumulative / k_values
    spearman_reason = None
    for k in range(K_MIN_FOR_SPEARMAN, k_max + 1):
        predicted_ranks = rankdata(ordered_prediction[:k])
        truth_ranks = rankdata(ordered_delta[:k])
        if np.ptp(predicted_ranks) == 0:
            spearman_reason = spearman_reason or CONSTANT_PREDICTION
            continue
        if np.ptp(truth_ranks) == 0:
            spearman_reason = spearman_reason or CONSTANT_TRUTH
            continue
        spearman[k - 1] = np.corrcoef(predicted_ranks, truth_ranks)[0, 1]
    if np.isnan(spearman).all():
        reasons['Spearman@k'] = spearman_reason or K_BELOW_THREE
    elif np.isnan(spearman[:K_MIN_FOR_SPEARMAN - 1]).any():
        reasons['Spearman@k'] = K_BELOW_THREE
    qualifying = np.flatnonzero(
        np.isfinite(spearman) & (precision >= PRECISION_TARGET) & (spearman >= SPEARMAN_TARGET))
    if len(qualifying):
        max_k = float(k_values[qualifying[-1]])
        recall = float(cumulative[int(max_k) - 1] / matched_total)
    else:
        reasons['max_k@P80S40'] = NO_QUALIFYING_K
        reasons['Recall@max_k'] = NO_QUALIFYING_K
    return GroupIndirectMetrics(k_values, precision, spearman, max_k, recall,
                                matched_total, reasons)


@dataclass(frozen=True)
class DirectAggregate:
    """Run-level direct metric mean over defined test wells, with counts and reasons."""
    values: dict
    counts: dict
    reasons: dict


def aggregate_direct(wells):
    """Equal-weight mean of per-well direct metrics over the defined wells."""
    wells = list(wells)
    values = {}
    counts = {'wells': len(wells)}
    reasons = {}
    for name in DIRECT_METRICS:
        samples = [well.values[name] for well in wells if np.isfinite(well.values[name])]
        counts[name] = len(samples)
        if samples:
            values[name] = float(np.mean(samples))
        else:
            values[name] = float('nan')
            reasons[name] = next((well.reasons[name] for well in wells if name in well.reasons),
                                 NO_DEFINED_UNITS)
    return DirectAggregate(values=values, counts=counts, reasons=reasons)


@dataclass(frozen=True)
class IndirectAggregate:
    """Run-level indirect metric means over defined test groups, with counts and reasons."""
    k: np.ndarray
    precision: np.ndarray
    spearman: np.ndarray
    max_k: float
    recall: float
    counts: dict
    reasons: dict


def aggregate_indirect(groups):
    """Equal-weight macro mean over defined test groups; no well weighting."""
    groups = list(groups)
    if not groups:
        raise ValueError('at least one group result is required')
    k = np.asarray(groups[0].k)
    if any(not np.array_equal(np.asarray(group.k), k) for group in groups):
        raise ValueError('group results must share one k axis')
    counts = {'groups': len(groups)}
    reasons = {}

    def curve(attribute, name):
        matrix = np.vstack([np.asarray(getattr(group, attribute), dtype=np.float64)
                            for group in groups])
        defined = np.isfinite(matrix)
        count = defined.sum(axis=0)
        mean = np.full(matrix.shape[1], np.nan)
        np.divide(np.where(defined, matrix, 0.).sum(axis=0), count, out=mean, where=count > 0)
        counts[name] = count
        if not count.any():
            reasons[name] = NO_DEFINED_GROUPS
        return mean

    def scalar(attribute, name):
        samples = [value for group in groups
                   if np.isfinite(value := float(getattr(group, attribute)))]
        counts[name] = len(samples)
        if not samples:
            reasons[name] = NO_DEFINED_GROUPS
            return float('nan')
        return float(np.mean(samples))

    precision = curve('precision', 'Precision@k')
    spearman = curve('spearman', 'Spearman@k')
    max_k = scalar('max_k', 'max_k@P80S40')
    recall = scalar('recall', 'Recall@max_k')
    return IndirectAggregate(k=k, precision=precision, spearman=spearman, max_k=max_k,
                             recall=recall, counts=counts, reasons=reasons)
