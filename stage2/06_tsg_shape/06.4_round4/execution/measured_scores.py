"""Pair held-out scores by frozen fold and training seed."""

from itertools import product

import numpy as np
import pandas as pd

PRIMARY_METRICS = ("Δ PCC", "Δ RMSE")
CONTRAST_METRICS = (*PRIMARY_METRICS, "Δ Spearman", "Δ E-distance")
CONTRASTS = {
    "XΔ-XΔΔ − XX-XXΔ": ("XΔ-XΔΔ", "XX-XXΔ"),
    "XΔ-XΔ − XΔ-XΔΔ": ("XΔ-XΔ", "XΔ-XΔΔ"),
}
SCORE_FIELDS = ["shape", "fold", "seed", "metric"]


def paired_report(scores, coordinates, metrics):
    """Require one value or an explicit undefined score for every metric and coordinate."""
    expected = pd.MultiIndex.from_tuples(
        [(*coordinate, metric) for coordinate, metric in product(coordinates, metrics)],
        names=SCORE_FIELDS,
    )
    indexed = scores.set_index(SCORE_FIELDS)
    if (not expected.is_unique or not indexed.index.is_unique or len(indexed) != len(expected)
            or not expected.difference(indexed.index).empty):
        raise ValueError("score rows differ from the registered arm, fold, seed and metric roster")
    counts = indexed[["n_conditions", "n_eligible", "n_undefined"]]
    if (counts.isna().any().any() or (counts < 0).any().any()
            or not counts.eq(counts.round()).all().all()
            or not counts.n_conditions.eq(counts.n_eligible + counts.n_undefined).all()
            or not indexed.value.notna().eq(counts.n_eligible.gt(0)).all()
            or not np.isfinite(indexed.value.dropna()).all()):
        raise ValueError("invalid score eligibility or nonfinite defined value")

    paired = indexed.reindex(expected).reset_index()
    pair_columns = ["fold", "seed", "metric"]
    by_arm = paired.pivot(index=pair_columns, columns="shape", values="value")
    paired["control_score"] = by_arm["XX-XX"].reindex(pd.MultiIndex.from_frame(
        paired[pair_columns])).to_numpy()
    paired["score_difference"] = paired.value - paired.control_score
    summary = paired.groupby(["shape", "metric"], sort=False, as_index=False).agg(
        value=("value", "mean"),
        sample_sd=("value", "std"),
        n_defined=("value", "count"),
        n_coordinates=("value", "size"),
        n_eligible_conditions=("n_eligible", "sum"),
        n_undefined_conditions=("n_undefined", "sum"),
        score_difference=("score_difference", "mean"),
        score_difference_sd=("score_difference", "std"),
        n_paired=("score_difference", "count"),
    )

    positions = tuple(dict.fromkeys((fold, seed) for _, fold, seed in coordinates))
    contrast_index = [(fold, seed, metric)
                      for (fold, seed), metric in product(positions, CONTRAST_METRICS)]
    focus = by_arm.reindex(contrast_index)
    contrasts = pd.concat([
        pd.DataFrame({
            "left_value": focus[left], "right_value": focus[right],
            "value": focus[left] - focus[right],
        }).reset_index().assign(contrast=label, left_shape=left, right_shape=right)
        for label, (left, right) in CONTRASTS.items()
    ], ignore_index=True)
    aggregates = {
        "value": ("value", "mean"),
        "sample_sd": ("value", "std"),
        "n_defined": ("value", "count"),
        "n_coordinates": ("value", "size"),
    }
    by_fold = contrasts.groupby(["contrast", "metric", "fold"], sort=False, as_index=False).agg(**aggregates)
    overall = contrasts.groupby(["contrast", "metric"], sort=False, as_index=False).agg(**aggregates)
    return paired, summary, contrasts, by_fold, overall
