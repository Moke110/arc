"""Pair held-out scores and summarize the registered Round 3 contrasts."""

from itertools import product

import numpy as np
import pandas as pd

PRIMARY_METRICS = ("Δ PCC", "Δ RMSE")
CONTRASTS = {
    "4−1": (4, 1),
    "4−3": (4, 3),
    "5−2": (5, 2),
    "5−3": (5, 3),
    "4−1−3+0": (4, 1, 3, 0),
    "5−2−3+0": (5, 2, 3, 0),
}
COORDINATE_FIELDS = ["shape", "fold", "seed"]
SCORE_FIELDS = [*COORDINATE_FIELDS, "metric"]


def paired_report(scores, coordinates, metrics):
    """Require the complete roster and pair only observations from the same fold and seed."""
    expected = pd.MultiIndex.from_tuples(
        [(*coordinate, metric) for coordinate, metric in product(coordinates, metrics)],
        names=SCORE_FIELDS,
    )
    indexed = scores.set_index(SCORE_FIELDS)
    if (not indexed.index.is_unique or len(indexed) != len(expected)
            or not expected.difference(indexed.index).empty):
        raise ValueError("score rows differ from the registered arm, fold, seed and metric roster")
    counts = indexed[["n_conditions", "n_eligible", "n_undefined"]]
    if (counts.isna().any().any() or (counts < 0).any().any()
            or not counts.n_conditions.eq(counts.n_eligible + counts.n_undefined).all()
            or not indexed.value.notna().eq(counts.n_eligible.gt(0)).all()
            or not np.isfinite(indexed.value.dropna()).all()):
        raise ValueError("invalid score eligibility or nonfinite defined value")
    paired = indexed.reindex(expected).reset_index()
    shapes = list(dict.fromkeys(shape for shape, _, _ in coordinates))
    coordinate_keys = ["fold", "seed", "metric"]
    scores_by_arm = paired.pivot(index=coordinate_keys, columns="shape", values="value")
    paired["control_score"] = scores_by_arm[shapes[0]].reindex(pd.MultiIndex.from_frame(
        paired[coordinate_keys])).to_numpy()
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
    positions = list(dict.fromkeys((fold, seed) for _, fold, seed in coordinates))
    primary_index = [(fold, seed, metric) for (fold, seed), metric in product(positions, PRIMARY_METRICS)]
    primary = scores_by_arm.reindex(primary_index)
    contrast_tables = []
    for label, (combined, first, *other) in CONTRASTS.items():
        difference = primary[shapes[combined]] - primary[shapes[first]]
        if other:
            second, original = other
            difference = difference - primary[shapes[second]] + primary[shapes[original]]
        contrast_tables.append(difference.rename("value").reset_index().assign(contrast=label))
    contrasts = pd.concat(contrast_tables, ignore_index=True)
    contrast_summary = contrasts.groupby(["contrast", "metric"], sort=False, as_index=False).agg(
        value=("value", "mean"),
        sample_sd=("value", "std"),
        n_defined=("value", "count"),
        n_coordinates=("value", "size"),
    )
    return paired, summary, contrasts, contrast_summary
