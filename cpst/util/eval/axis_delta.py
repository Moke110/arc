"""Held-out signed Delta scores on the full prediction axis."""

import numpy as np
import pandas as pd

from cpst_util.metrics import signed_jaccard_scores
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
