"""Signed condition scores and empirical fold-level Delta distribution metrics."""

import math
from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist

from cpst_util.metrics import vector_metrics

METRICS = ("JAC_top200", "JAC_pos100", "JAC_neg100", "Δ PCC", "Δ Spearman",
           "Δ RMSE", "Δ E-distance", "Δ MMD", "Δ R²")
DISTRIBUTION = METRICS[6:8]
CONDITION = (*METRICS[:6], METRICS[-1])


def _extreme(values, genes, sign, limit):
    magnitude = sign * values
    nonzero = np.flatnonzero(magnitude > 0)
    order = np.lexsort((genes[nonzero], -magnitude[nonzero]))
    return set(genes[nonzero[order[:limit]]])


def _jaccard(hits, total):
    return hits / total if total else math.nan


def signed_jaccard(truth, predicted, gene_indices, *, limit=100):
    """Return signed union, positive and negative Jaccard, breaking ties by gene index."""
    truth, predicted, genes = map(np.asarray, (truth, predicted, gene_indices))
    if truth.ndim != 1 or truth.shape != predicted.shape or genes.shape != truth.shape:
        raise ValueError("Jaccard requires paired vectors on one ordered gene axis")
    hits, totals = [], []
    for sign in (1, -1):
        observed = _extreme(truth, genes, sign, limit)
        proposed = _extreme(predicted, genes, sign, limit)
        hits.append(len(observed & proposed))
        totals.append(len(observed | proposed))
    return (_jaccard(sum(hits), sum(totals)), *map(_jaccard, hits, totals))


def condition_scores(truth, predicted, gene_indices):
    """Score one condition with the existing Stage 2 vector-metric semantics."""
    vector = vector_metrics(truth, predicted)
    values = (*signed_jaccard(truth, predicted, gene_indices),
              *(vector[key] for key in ("pcc", "spearman", "rmse", "r2")))
    return dict(zip(CONDITION, values, strict=True))


class DiseaseScore(NamedTuple):
    f1: float
    tp: int
    fp: int
    fn: int


def signed_disease_f1(predicted, gene_indices, cm_direction, disease_genes):
    """Compare all signed predictions against nonzero CM directions among the disease genes."""
    prediction, genes = np.asarray(predicted), np.asarray(gene_indices)
    observed, disease = np.asarray(cm_direction), np.asarray(disease_genes)
    prediction_sign = np.sign(prediction)
    present = prediction_sign != 0
    common = observed[genes]
    matches = present & (prediction_sign == common) & np.isin(genes, disease)
    tp = int(np.count_nonzero(matches))
    fp = int(np.count_nonzero(present)) - tp
    fn = int(np.count_nonzero(observed[disease])) - tp
    return DiseaseScore(2 * tp / (2 * tp + fp + fn) if tp + fp + fn else math.nan, tp, fp, fn)


def fold_bandwidth(training_delta):
    """Median scaled Euclidean distance of unordered distinct HVG-DEG 5k training row pairs."""
    values = np.asarray(training_delta)
    if values.ndim != 2 or values.shape[0] < 2 or not values.shape[1] or not np.isfinite(values).all():
        raise ValueError("bandwidth requires at least two finite training Delta rows")
    bandwidth = float(np.median(pdist(values, "euclidean") / math.sqrt(values.shape[1])))
    if bandwidth <= 0:
        raise ValueError("median training pair distance must be positive for the RBF bandwidth")
    return bandwidth


def distribution_scores(truth, predicted, bandwidth):
    """Return scaled E-distance and square root of biased empirical RBF MMD²."""
    truth, predicted = np.asarray(truth), np.asarray(predicted)
    if (truth.ndim != 2 or not truth.size or truth.shape != predicted.shape
            or not np.isfinite(truth).all() or not np.isfinite(predicted).all()
            or not math.isfinite(bandwidth) or bandwidth <= 0):
        raise ValueError("distribution metrics require paired finite condition matrices and a positive bandwidth")

    def means(left, right):
        distances = cdist(left, right, "euclidean") / math.sqrt(truth.shape[1])
        kernel = np.exp(-0.5 * np.square(distances / bandwidth))
        return math.fsum(distances.flat) / distances.size, math.fsum(kernel.flat) / kernel.size

    cross_d, cross_k = means(truth, predicted)
    truth_d, truth_k = means(truth, truth)
    predicted_d, predicted_k = means(predicted, predicted)
    return (max(0.0, float(2 * cross_d - (truth_d + predicted_d))),
            math.sqrt(max(0.0, float(truth_k + predicted_k - 2 * cross_k))))


def summarize_conditions(scores):
    """Average defined condition scores and keep totals for each eligibility class."""
    scores = np.asarray(scores)
    defined = scores[np.isfinite(scores)]
    return {"value": math.fsum(defined) / len(defined) if len(defined) else math.nan,
            "n_conditions": len(scores), "n_eligible": len(defined),
            "n_undefined": len(scores) - len(defined)}


def summarize_folds(records, keys):
    """Average defined folds equally while retaining condition and fold counts."""
    result = []
    counts = ["n_conditions", "n_eligible", "n_undefined"]
    for group, fold in records.groupby(keys, sort=False):
        scores = summarize_conditions(fold.value)
        result.append(dict(zip(keys, group, strict=True)) | fold[counts].sum().to_dict() | {
            "value": scores["value"], "n_folds": scores["n_conditions"],
            "n_eligible_folds": scores["n_eligible"], "n_undefined_folds": scores["n_undefined"],
        })
    return pd.DataFrame(result, columns=[*keys, "value", "n_folds", "n_eligible_folds",
                                         "n_undefined_folds", *counts])
