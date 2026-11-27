"""Generic expression-vector metrics shared by Stage 2 studies."""

import math

import numpy as np
from scipy.stats import rankdata


def vector_metrics(truth, prediction):
    """Score paired vectors with the accepted constant-truth mask."""
    truth = np.asarray(truth, dtype=np.float64)
    prediction = np.asarray(prediction, dtype=np.float64)
    if (truth.ndim != 1 or truth.size == 0 or truth.shape != prediction.shape
            or not np.isfinite(truth).all() or not np.isfinite(prediction).all()):
        raise ValueError("expected paired nonempty, finite vectors")

    constant_truth = bool(np.ptp(truth) == 0)
    constant_prediction = bool(np.ptp(prediction) == 0)
    squared_error = (truth - prediction) ** 2
    denominator = float(np.sum((truth - truth.mean()) ** 2))
    if constant_truth:
        pcc = spearman = float("nan")
    elif constant_prediction:
        pcc = spearman = 0.0
    else:
        pcc = float(np.corrcoef(truth, prediction)[0, 1])
        spearman = float(np.corrcoef(rankdata(truth), rankdata(prediction))[0, 1])
    return {
        "pcc": pcc,
        "spearman": spearman,
        "rmse": math.sqrt(float(squared_error.mean())),
        "r2": 1 - float(squared_error.sum()) / denominator if denominator else float("nan"),
        "constant_truth": constant_truth,
        "constant_prediction": constant_prediction,
    }


JACCARD_CUTOFFS = (200, 500, 1000, 2000)
JACCARD_WEIGHTS = (0.125, 0.25, 0.5, 1.0)


def _rank_by_sign(values, genes, sign):
    magnitude = sign * values
    nonzero = np.flatnonzero(magnitude > 0)
    order = np.lexsort((genes[nonzero], -magnitude[nonzero]))
    return genes[nonzero[order]]


def _signed_rankings(truth, predicted, gene_indices):
    truth, predicted, genes = map(np.asarray, (truth, predicted, gene_indices))
    if truth.ndim != 1 or truth.shape != predicted.shape or genes.shape != truth.shape:
        raise ValueError("Jaccard requires paired vectors on one ordered gene axis")
    return [(_rank_by_sign(truth, genes, sign), _rank_by_sign(predicted, genes, sign))
            for sign in (1, -1)]


def _jaccard(hits, total):
    return hits / total if total else math.nan


def _direction_counts(rankings, limit):
    for truth, predicted in rankings:
        observed, proposed = set(truth[:limit]), set(predicted[:limit])
        yield len(observed & proposed), len(observed | proposed)


def signed_jaccard(truth, predicted, gene_indices, *, limit=100):
    """Return combined, positive and negative Jaccard at one directional limit."""
    rankings = _signed_rankings(truth, predicted, gene_indices)
    hits, totals = zip(*_direction_counts(rankings, limit), strict=True)
    return _jaccard(sum(hits), sum(totals)), *map(_jaccard, hits, totals)


def signed_jaccard_scores(truth, predicted, gene_indices):
    """Return JAC-top200 and wJAC-2000 for one matched-control Delta pair."""
    rankings = _signed_rankings(truth, predicted, gene_indices)
    scores = []
    for cutoff in JACCARD_CUTOFFS:
        hits, totals = zip(*_direction_counts(rankings, cutoff // 2), strict=True)
        scores.append(_jaccard(sum(hits), sum(totals)))
    weighted = math.fsum(weight * score for weight, score in zip(JACCARD_WEIGHTS, scores, strict=True))
    return scores[0], weighted / math.fsum(JACCARD_WEIGHTS)
