"""Held-out signed Delta scores on the full prediction axis."""

import numpy as np
import pandas as pd

from axis_protocol import REFERENCE_AXIS
from stage2_util.metrics import signed_jaccard_scores
from tsg_axis.metrics import CONDITION, DISTRIBUTION, METRICS, condition_scores, distribution_scores, summarize_conditions

WEIGHTED = "wJAC-2000"
REPORT_METRICS = (*METRICS, WEIGHTED)
REPORT_CONDITION = (*CONDITION, WEIGHTED)


def score_deltas(measured, predicted, genes, bandwidth):
    """Score each held-out condition on its registered prediction axis."""
    if (measured.ndim != 2 or measured.shape != predicted.shape
            or measured.shape[1] != len(genes) or len(set(genes)) != len(genes)):
        raise ValueError("scoring requires aligned Delta matrices and unique ordered genes")
    records = []
    for actual, estimate in zip(measured, predicted, strict=True):
        scores = condition_scores(actual, estimate, genes)
        top200, weighted = signed_jaccard_scores(actual, estimate, genes)
        if scores["JAC_top200"] != top200 and not (np.isnan(scores["JAC_top200"]) and np.isnan(top200)):
            raise ValueError("full-axis JAC-top200 differs from the accepted scoring kernel")
        records.append({**scores, WEIGHTED: weighted})
    conditions = pd.DataFrame(records)
    distribution = dict(zip(DISTRIBUTION, distribution_scores(measured, predicted, bandwidth), strict=True))
    by_fold = []
    for metric in REPORT_METRICS:
        summary = summarize_conditions(conditions[metric]) if metric in REPORT_CONDITION else {
            "value": distribution[metric], "n_conditions": len(conditions),
            "n_eligible": len(conditions), "n_undefined": 0,
        }
        by_fold.append({"metric": metric, "width": len(genes), **summary})
    return conditions, pd.DataFrame(by_fold)


def paired_scores(scores, roster, *, metrics=REPORT_METRICS):
    """Pair coordinates by fold and seed and summarize equal-weight defined scores."""
    keys = ["axis", "fold", "seed", "metric"]
    expected = pd.MultiIndex.from_tuples(
        [(*coordinate, metric) for coordinate in roster for metric in metrics], names=keys)
    indexed = scores.set_index(keys)
    if not indexed.index.is_unique or len(indexed) != len(expected) or not expected.difference(indexed.index).empty:
        raise ValueError("each coordinate must have each registered metric exactly once")
    counts = indexed[["n_conditions", "n_eligible", "n_undefined"]]
    if (counts.isna().any().any() or (counts < 0).any().any()
            or not counts.eq(counts.round()).all().all()
            or not counts.n_conditions.eq(counts.n_eligible + counts.n_undefined).all()
            or not indexed.value.notna().eq(counts.n_eligible.gt(0)).all()
            or not np.isfinite(indexed.value.dropna()).all()):
        raise ValueError("values must match the registered eligibility counts")
    paired = indexed.reindex(expected).reset_index()
    columns = ["fold", "seed", "metric"]
    reference = paired[paired.axis.eq(REFERENCE_AXIS)].set_index(columns)
    pairs = pd.MultiIndex.from_frame(paired[columns])
    if not pairs.difference(reference.index).empty:
        raise ValueError("each fold and seed needs a matching reference coordinate")
    paired["reference_value"] = reference.value.reindex(pairs).to_numpy()
    paired["paired_difference"] = paired.value - paired.reference_value
    summary = paired.groupby(["axis", "metric"], sort=False, as_index=False).agg(
        mean=("value", "mean"), sample_sd=("value", "std"), n_defined=("value", "count"),
        paired_mean=("paired_difference", "mean"), paired_sd=("paired_difference", "std"),
        n_paired=("paired_difference", "count"), n_coordinates=("value", "size"),
        n_eligible=("n_eligible", "sum"), n_undefined=("n_undefined", "sum"),
        width=("width", "first"))
    return paired, summary
