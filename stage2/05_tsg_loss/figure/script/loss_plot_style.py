"""Shared CSV bindings and styling for the Step 05 report panels."""

from collections import defaultdict
import math
from statistics import mean

from matplotlib import pyplot as plt
import numpy as np

from stage2_util import palette, plotting

LOWER_IS_BETTER = frozenset(("Δ RMSE", "Δ E-distance", "Δ MMD"))
FOLD_SCORES = "fold-scores.csv"


def metric_direction(axis, metric, label):
    """Label the preferred direction of the metric."""
    lower = metric in LOWER_IS_BETTER
    axis.set_title(f"{label} {'↓' if lower else '↑'}", fontweight="bold")


def read_replicates(path, scores, *, coordinates, reference, folds):
    """Bind each reported group to its raw and paired fold/seed observations."""
    groups = defaultdict(dict)
    for row in plotting.read_rows(path):
        key = tuple(row[name] for name in (*coordinates, "metric"))
        run = int(row["fold"]), int(row["seed"])
        value = float(row["value"]) if row["value"] else math.nan
        if math.isinf(value) or run in groups[key]:
            raise ValueError(f"invalid fold score for {key}, run {run}")
        groups[key][run] = value
    if groups.keys() != scores.keys():
        raise ValueError("fold scores must contain the complete metric grid")

    first = next(iter(groups.values()))
    if {fold for fold, _ in first} != set(folds) or any(rows.keys() != first.keys()
                                                        for rows in groups.values()):
        raise ValueError("fold scores must contain the same fold/seed roster for every group")

    replicates = {}
    for key, rows in groups.items():
        reference_rows = groups[(*reference, key[-1])]
        values = [rows[run] for run in sorted(rows)]
        differences = [rows[run] - reference_rows[run] for run in sorted(rows)]
        defined = [value for value in values if math.isfinite(value)]
        average = mean(defined) if defined else math.nan
        expected = scores[key]["value"]
        if not (math.isclose(average, expected, rel_tol=1e-12, abs_tol=1e-12)
                or math.isnan(average) and math.isnan(expected)):
            raise ValueError(f"fold scores disagree with equal-fold summary for {key}")
        replicates[key] = {"value": values, "score_difference": differences}
    return replicates


def score_bar(axis, position, values, *, hue, width, depth=2):
    """Draw the replicate score group and label undefined groups."""
    bars, _, _ = plotting.replicate_bar(
        axis, position, values, width=width, hue=hue, depth=depth, linewidth=.8,
    )
    if bars is None:
        axis.text(position, .025, "NA", ha="center", fontsize=7,
                  transform=axis.get_xaxis_transform(), zorder=4)


def setup_score_axis(axis, x, labels, *, difference):
    axis.set(xticks=x, xticklabels=labels, xlim=(x[0] - .5, x[-1] + .5))
    axis.axhline(0, color="#666666", linewidth=.8, zorder=1.5)
    axis.set_axisbelow(True)
    axis.grid(axis="y", color="#777777", alpha=.25)
    axis.spines[["top", "right"]].set_visible(False)
    if difference:
        limit = max(map(abs, axis.get_ylim()))
        axis.set_ylim(-limit, limit)


def render_weight_panels(replicates, field, output, alphas, labels,
                         *, xlabel, reference, rotate=False):
    """Render one score panel per metric for a one-dimensional weight sweep."""
    difference = field == "score_difference"
    fig, axes = plt.subplots(3, 3, figsize=(13, 9), layout="constrained")
    x = np.arange(len(alphas))
    for (metric, label), axis in zip(labels.items(), axes.flat, strict=True):
        metric_direction(axis, metric, label)
        for position, alpha in enumerate(alphas):
            score_bar(axis, position, replicates[alpha, metric][field],
                      hue=palette.HUES[position], width=.68)
        setup_score_axis(axis, x, alphas, difference=difference)
        if rotate:
            for tick in axis.get_xticklabels():
                tick.set(rotation=45, ha="right", fontsize=9)

    fig.supxlabel(xlabel)
    fig.supylabel(f"Paired score difference vs. α = {reference}" if difference else "Equal-fold score")
    plotting.save_figure(fig, output, dpi=180)
    plt.close(fig)
