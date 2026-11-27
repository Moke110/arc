"""Held-out full-axis and disease-restricted signed Delta scores."""

import numpy as np
import pandas as pd

from axis_protocol import REFERENCE_AXIS, SCOPES
from tsg_axis.metrics import CONDITION, DISTRIBUTION, METRICS, condition_scores, distribution_scores, summarize_conditions

COVERAGE = "direction-correct coverage"


def score_deltas(measured, predicted, genes, union, denominator, bandwidths):
    """Score the two registered gene scopes and direction-correct disease coverage."""
    if (measured.shape != predicted.shape or measured.ndim != 2
            or measured.shape[1] != len(genes) or len(set(genes)) != len(genes)
            or denominator < len(union) or set(bandwidths) != set(SCOPES)):
        raise ValueError("score scopes require aligned Delta matrices, genes and fold bandwidths")
    positions = np.isin(genes, list(union))
    condition_tables, scores = [], []
    for scope, keep in zip(SCOPES, (slice(None), positions), strict=True):
        truth, proposal, scope_genes = measured[:, keep], predicted[:, keep], genes[keep]
        if not len(scope_genes):
            raise ValueError("empty prediction gene scope")
        conditions = pd.DataFrame(condition_scores(actual, estimate, scope_genes)
                                  for actual, estimate in zip(truth, proposal, strict=True))
        distribution = dict(zip(DISTRIBUTION, distribution_scores(truth, proposal, bandwidths[scope]), strict=True))
        metrics = METRICS
        if scope == "disease":
            signs = np.sign(truth)
            correct = np.count_nonzero((signs != 0) & (signs == np.sign(proposal)), axis=1)
            conditions[COVERAGE] = correct * (100 / denominator)
            metrics = (*METRICS, COVERAGE)
        for metric in metrics:
            summary = summarize_conditions(conditions[metric]) if metric in (*CONDITION, COVERAGE) else {
                "value": distribution[metric], "n_conditions": len(conditions),
                "n_eligible": len(conditions), "n_undefined": 0,
            }
            scores.append({"scope": scope, "metric": metric, "overlap": len(scope_genes), **summary})
        condition_tables.append(conditions.assign(scope=scope, held_out_index=np.arange(len(conditions))))
    return pd.concat(condition_tables, ignore_index=True), pd.DataFrame(scores)


def paired_scores(scores, roster, *, metrics=None):
    """Pair coordinates by fold and seed and summarize equal-weight defined scores."""
    if metrics is None:
        metrics = {"axis": METRICS, "disease": (*METRICS, COVERAGE)}
    keys = ["axis", "fold", "seed", "scope", "metric"]
    expected = pd.MultiIndex.from_tuples(
        [(*coordinate, scope, metric) for coordinate in roster
         for scope, names in metrics.items() for metric in names], names=keys)
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
    columns = ["fold", "seed", "scope", "metric"]
    reference = paired[paired.axis.eq(REFERENCE_AXIS)].set_index(columns)
    pairs = pd.MultiIndex.from_frame(paired[columns])
    if not pairs.difference(reference.index).empty:
        raise ValueError("each fold and seed needs a matching reference coordinate")
    paired["reference_value"] = reference.value.reindex(pairs).to_numpy()
    paired["paired_difference"] = paired.value - paired.reference_value
    summary = paired.groupby(["axis", "scope", "metric"], sort=False, as_index=False).agg(
        mean=("value", "mean"), sample_sd=("value", "std"), n_defined=("value", "count"),
        paired_mean=("paired_difference", "mean"), paired_sd=("paired_difference", "std"),
        n_paired=("paired_difference", "count"), n_coordinates=("value", "size"),
        n_eligible=("n_eligible", "sum"), n_undefined=("n_undefined", "sum"),
        overlap=("overlap", "first"))
    return paired, summary
