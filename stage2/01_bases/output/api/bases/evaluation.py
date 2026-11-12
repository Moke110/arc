"""Stage 2 evaluation semantics shared by every model-comparison step.

Owns the measured-effect gate, the per-well direct metrics, the predicted
absolute rank table, the group-level ``@d`` reliability family and run-level
aggregation. Every function is pure over arrays plus identity metadata: no file
I/O, no cached state.

The reliability family is total by convention (ADR 0064): the four
empty/non-empty window/measured-effect cells take declared values, so no
reliability metric is NaN and none carries a reason or a defined-group count.
The retired ``@k`` family, its ``|log2FC| >= 2`` gate and its NaN-plus-reason
machinery are gone. The retired AUPRC/PR implementation in ``bases.metrics`` is
deliberately not imported or modified here.
"""
from dataclasses import dataclass
import numpy as np

GENE_COUNT = 18164
DIRECT_METRICS = ('PCC', 'Spearman', 'RMSE', 'R2')
RELIABILITY_METRICS = ('Precision@d', 'Recall@d', 'F1@d', 'Spearman@d', 'Reliability@d')

PADJ_THRESHOLD = 0.05
THRESHOLD_GRID_POINTS = 100
MIN_HITS_FOR_SPEARMAN = 3
SPEARMAN_NEUTRAL = 1.0

# Undefined direct-metric reasons; the reliability family never uses them.
CONSTANT_PREDICTION = 'constant_prediction'
CONSTANT_TRUTH = 'constant_truth'
NO_DEFINED_UNITS = 'no_defined_units'


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
class MeasuredEffects:
    """Full-axis measured-effect gate: direction from ``sign(log2FC)``, 0 elsewhere.

    ``magnitude`` is ``|log2FC|`` on the measured genes and 0 elsewhere, so the
    ordering term always has a finite measured magnitude to rank.
    """
    signed: np.ndarray
    magnitude: np.ndarray
    indices: np.ndarray

    @property
    def count(self):
        return int(len(self.indices))


def measured_effect_gate(padj, log2_fold_change):
    """Measured effects: ``padj < 0.05`` with direction ``sign(log2FC)`` (ADR 0064).

    There is no fold-change magnitude gate. A NaN ``padj`` is not significant; a
    selected gene must have a finite, nonzero ``log2FC``.
    """
    arrays = _shared_axis(padj=padj, log2_fold_change=log2_fold_change)
    padj = arrays['padj'].astype(np.float64)
    effect = arrays['log2_fold_change'].astype(np.float64)
    finite = np.isfinite(padj)
    if np.any(np.isinf(padj)) or np.any((padj[finite] < 0) | (padj[finite] > 1)):
        raise ValueError('invalid adjusted p-values')
    selected = padj < PADJ_THRESHOLD
    if np.any(selected & (~np.isfinite(effect) | (effect == 0))):
        raise ValueError('significant gene has invalid direction')
    signed = np.zeros(len(padj), dtype=np.int8)
    magnitude = np.zeros(len(padj), dtype=np.float64)
    indices = np.flatnonzero(selected).astype(np.int64)
    signed[indices] = np.sign(effect[indices]).astype(np.int8)
    magnitude[indices] = np.abs(effect[indices])
    return MeasuredEffects(signed=signed, magnitude=magnitude, indices=indices)


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
class ThresholdGrid:
    """One run's declared 100-point threshold grid.

    ``d_prime = k/100`` is the scale-free plotted axis; ``thresholds = d'·M`` is
    the predicted-magnitude threshold in the run's own units, where ``M`` is the
    run's largest ``|Δ̂_mean|`` over its test groups and the whole gene axis.
    """
    d_prime: np.ndarray
    thresholds: np.ndarray
    max_absolute_delta: float


def threshold_grid(max_absolute_delta):
    """Declared scan ``d = k·0.01·M`` for ``k = 1..100`` (ADR 0064), as ``d'·M``.

    ``M`` must be a positive finite magnitude, so ``d' = d/M = k/100`` and the
    strictest point is exactly ``M``. Window membership is ``|Δ̂| >= d``.
    """
    normaliser = float(max_absolute_delta)
    if not np.isfinite(normaliser) or normaliser <= 0:
        raise ValueError('max_absolute_delta must be a positive finite magnitude')
    steps = np.arange(1, THRESHOLD_GRID_POINTS + 1, dtype=np.float64)
    d_prime = steps / THRESHOLD_GRID_POINTS
    return ThresholdGrid(d_prime=d_prime, thresholds=d_prime * normaliser,
                         max_absolute_delta=normaliser)


def precision_at_d(n_window, n_hits, n_measured):
    """``|H|/|W|`` with the declared empty-window totals.

    An empty window scores 1 when the group has no measured effect and 0 when it
    has one; a non-empty window scores ``|H|/|W|``.
    """
    if n_window == 0:
        return 1.0 if n_measured == 0 else 0.0
    return n_hits / n_window


def recall_at_d(n_measured, n_hits, n_window):
    """``|H|/|measured effects|`` with the declared empty-label totals.

    A group with no measured effect scores 1 when its window is also empty and 0
    when it emitted a window; otherwise the recall is ``|H|/|measured|``.
    """
    if n_measured == 0:
        return 1.0 if n_window == 0 else 0.0
    return n_hits / n_measured


def f1_at_d(precision, recall):
    """Harmonic mean of ``Precision@d`` and ``Recall@d``; ``0`` when both are 0.

    Harmonic averaging is load-bearing: an arithmetic mean peaks at the smallest
    window and would make the reliability peak select a one-gene table.
    """
    total = precision + recall
    if total == 0:
        return 0.0
    return 2.0 * precision * recall / total


def spearman_at_d(predicted_magnitude, measured_magnitude):
    """Spearman inside the direction-correct hits: predicted vs measured magnitude.

    Neutral (``1.0``, so the reliability multiplier is 1) below
    ``MIN_HITS_FOR_SPEARMAN`` hits or when either hit magnitude has no rank
    variation, because the ordering cannot be estimated there. Ties use average
    ranks.
    """
    predicted = np.asarray(predicted_magnitude, dtype=np.float64)
    measured = np.asarray(measured_magnitude, dtype=np.float64)
    if predicted.ndim != 1 or measured.ndim != 1 or len(predicted) != len(measured):
        raise ValueError('spearman magnitudes must share one hit axis')
    if not np.isfinite(predicted).all() or not np.isfinite(measured).all():
        raise ValueError('spearman magnitudes must be finite')
    if len(predicted) < MIN_HITS_FOR_SPEARMAN:
        return SPEARMAN_NEUTRAL
    from scipy.stats import rankdata

    predicted_ranks = rankdata(predicted)
    measured_ranks = rankdata(measured)
    if np.ptp(predicted_ranks) == 0 or np.ptp(measured_ranks) == 0:
        return SPEARMAN_NEUTRAL
    return float(np.corrcoef(predicted_ranks, measured_ranks)[0, 1])


def reliability_at_d(f1, spearman):
    """``F1@d · (1 + Spearman@d)/2``, in ``[0, 1]``."""
    return f1 * (1.0 + spearman) / 2.0


@dataclass(frozen=True)
class GroupReliabilityCurve:
    """One group's ``@d`` curve over a run's threshold grid.

    ``spearman`` is neutral (1.0) where the ordering term cannot be estimated, so
    every array is finite and ``reliability = f1 · (1 + spearman)/2`` holds at
    every point. ``n_window`` and ``n_hits`` vary per threshold; ``n_measured``
    is the group's fixed measured-effect count.
    """
    grid: ThresholdGrid
    precision: np.ndarray
    recall: np.ndarray
    f1: np.ndarray
    spearman: np.ndarray
    reliability: np.ndarray
    n_window: np.ndarray
    n_hits: np.ndarray
    n_measured: int

    @property
    def d_prime(self):
        return self.grid.d_prime

    @property
    def thresholds(self):
        return self.grid.thresholds


def reliability_curve(delta_mean, effects, grid):
    """Score one group's predicted mean Δ against its measured effects.

    The window is ``W(d) = {i : |Δ̂_mean| >= d}``; a direction-correct hit is a
    window gene that is a measured effect with matching direction. ``grid`` is
    the run-level grid, so a weak group's window can empty at high thresholds.
    """
    delta = _float_vector(delta_mean, 'delta_mean')
    signed = np.asarray(effects.signed)
    magnitude = np.asarray(effects.magnitude, dtype=np.float64)
    if signed.ndim != 1 or len(signed) != len(delta) or len(magnitude) != len(delta):
        raise ValueError('reliability curve requires one shared gene axis')
    thresholds = np.asarray(grid.thresholds, dtype=np.float64)
    if thresholds.ndim != 1 or len(thresholds) != THRESHOLD_GRID_POINTS:
        raise ValueError('threshold grid must hold 100 points')

    absolute = np.abs(delta)
    direction = np.sign(delta).astype(np.int8)
    measured = signed != 0
    n_measured = int(np.count_nonzero(measured))
    window = absolute[None, :] >= thresholds[:, None]
    hits = window & measured[None, :] & (direction[None, :] == signed[None, :])

    count = len(thresholds)
    precision = np.empty(count)
    recall = np.empty(count)
    f1 = np.empty(count)
    spearman = np.empty(count)
    reliability = np.empty(count)
    n_window = np.empty(count, dtype=np.int64)
    n_hits = np.empty(count, dtype=np.int64)
    for index in range(count):
        n_window[index] = int(np.count_nonzero(window[index]))
        n_hits[index] = int(np.count_nonzero(hits[index]))
        precision[index] = precision_at_d(int(n_window[index]), int(n_hits[index]), n_measured)
        recall[index] = recall_at_d(n_measured, int(n_hits[index]), int(n_window[index]))
        f1[index] = f1_at_d(float(precision[index]), float(recall[index]))
        hit_index = np.flatnonzero(hits[index])
        spearman[index] = spearman_at_d(absolute[hit_index], magnitude[hit_index])
        reliability[index] = reliability_at_d(float(f1[index]), float(spearman[index]))
    return GroupReliabilityCurve(
        grid=grid, precision=precision, recall=recall, f1=f1, spearman=spearman,
        reliability=reliability, n_window=n_window, n_hits=n_hits, n_measured=n_measured)


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
class ReliabilityAggregate:
    """Run-level ``@d`` macro curves over the evaluated fold's test groups.

    Every group contributes to every macro because the family is total, so
    ``counts['groups']`` is the population behind each metric curve; the
    ``null``/``changed`` split macros are means over the groups of each class.
    """
    grid: ThresholdGrid
    precision: np.ndarray
    recall: np.ndarray
    f1: np.ndarray
    spearman: np.ndarray
    reliability: np.ndarray
    null_reliability: np.ndarray
    changed_reliability: np.ndarray
    counts: dict

    @property
    def d_prime(self):
        return self.grid.d_prime

    @property
    def thresholds(self):
        return self.grid.thresholds


def aggregate_reliability(curves):
    """Equal-weight macro mean over one fold's test groups, plus the null split.

    ``curves`` must share one run's threshold grid. A class with no groups in
    that fold gets a zero curve and a zero count, never NaN.
    """
    curves = list(curves)
    if not curves:
        raise ValueError('at least one group curve is required')
    grid = curves[0].grid
    for curve in curves[1:]:
        if (not np.array_equal(curve.grid.d_prime, grid.d_prime)
                or not np.array_equal(curve.grid.thresholds, grid.thresholds)):
            raise ValueError('group curves must share one threshold grid')
    points = len(grid.thresholds)

    def macro(attribute):
        return np.mean(np.vstack([np.asarray(getattr(curve, attribute), dtype=np.float64)
                                  for curve in curves]), axis=0)

    precision, recall = macro('precision'), macro('recall')
    f1, spearman = macro('f1'), macro('spearman')
    reliability = macro('reliability')
    reliability_matrix = np.vstack([np.asarray(curve.reliability, dtype=np.float64)
                                    for curve in curves])
    measured = np.array([curve.n_measured for curve in curves], dtype=np.int64)
    null = measured == 0
    changed = ~null

    def split(mask):
        if not mask.any():
            return np.zeros(points)
        return reliability_matrix[mask].mean(axis=0)

    return ReliabilityAggregate(
        grid=grid, precision=precision, recall=recall, f1=f1, spearman=spearman,
        reliability=reliability, null_reliability=split(null),
        changed_reliability=split(changed),
        counts={'groups': len(curves), 'null_groups': int(null.sum()),
                'changed_groups': int(changed.sum())})
