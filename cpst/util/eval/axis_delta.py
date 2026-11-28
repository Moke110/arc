"""Held-out signed Delta scores on the full prediction axis."""

import pandas as pd

from tsg_axis.metrics import CONDITION, DISTRIBUTION, METRICS, condition_scores, distribution_scores, summarize_conditions

REPORT_METRICS = METRICS
REPORT_CONDITION = CONDITION


def score_deltas(measured, predicted, genes, bandwidth):
    """Score each held-out condition on its registered prediction axis."""
    if (measured.ndim != 2 or measured.shape != predicted.shape
            or measured.shape[1] != len(genes) or len(set(genes)) != len(genes)):
        raise ValueError("scoring requires aligned Delta matrices and unique ordered genes")
    conditions = pd.DataFrame(
        condition_scores(actual, estimate, genes)
        for actual, estimate in zip(measured, predicted, strict=True)
    )
    distribution = dict(zip(DISTRIBUTION, distribution_scores(measured, predicted, bandwidth), strict=True))
    by_fold = []
    for metric in REPORT_METRICS:
        summary = summarize_conditions(conditions[metric]) if metric in REPORT_CONDITION else {
            "value": distribution[metric], "n_conditions": len(conditions),
            "n_eligible": len(conditions), "n_undefined": 0,
        }
        by_fold.append({"metric": metric, "width": len(genes), **summary})
    return conditions, pd.DataFrame(by_fold)
