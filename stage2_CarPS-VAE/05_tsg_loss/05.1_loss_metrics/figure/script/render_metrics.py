"""Render equal-fold Delta scores and signed differences against native MSE."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

LOSS = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(LOSS.parent / "util/output/api"), str(LOSS / "figure/script")]

from loss_plot_style import (FOCUSED_Y_NOTE, FOLD_SCORES, metric_direction, read_replicates,
                             score_bar, setup_score_axis)
from stage2_util import palette, plotting

LOSSES = {"mse": "MSE", "cos": "Cos", "pcc": "PCC", "l1": "L1"}
ALPHAS = ("0.75", "1", "1.25")
LOSS_HUES = (220, 20, 140, 300)
ALPHA_DEPTHS = (0, 1, 3)
LABELS = {
    "JAC_top200": "JAC top200",
    "JAC_pos100": "JAC pos100",
    "JAC_neg100": "JAC neg100",
    "Δ PCC": "Δ PCC",
    "Δ Spearman": "Δ Spearman",
    "Δ RMSE": "Δ RMSE",
    "Δ E-distance": "Δ E-distance",
    "Δ MMD": "Δ MMD",
    "Δ R²": "Δ R²",
}
FIGURES = {"value": "delta-metrics.png",
           "score_difference": "delta-metrics-improvement.png"}


def read_scores(path):
    scores = {}
    for row in plotting.read_rows(path):
        key = row["loss"], row["alpha"], row["metric"]
        if key in scores:
            raise ValueError(f"duplicate score in complete loss grid: {key}")
        record = {name: float(row[name]) if row[name] else math.nan
                  for name in (*FIGURES, "native_mse_score")}
        if any(math.isinf(value) for value in record.values()):
            raise ValueError(f"infinite score for {key}")
        scores[key] = record
    if set(scores) != set(product(LOSSES, ALPHAS, LABELS)):
        raise ValueError("metric summary must contain the complete loss grid")
    for key, record in scores.items():
        baseline = scores["mse", "1", key[2]]["value"]
        actual = record["native_mse_score"], record["score_difference"]
        expected = baseline, record["value"] - baseline
        if any(a != b and not (math.isnan(a) and math.isnan(b))
               for a, b in zip(actual, expected, strict=True)):
            raise ValueError(f"signed score difference disagrees with the native MSE score: {key}")
    return scores


def render_panels(replicates, field, output):
    difference = field == "score_difference"
    fig, plots = plt.subplots(3, 3, figsize=(13, 9), layout="constrained")
    x = np.arange(len(LOSSES))
    width = .24
    for (metric, label), plot in zip(LABELS.items(), plots.flat, strict=True):
        metric_direction(plot, metric, label)
        groups = []
        for position, (form, hue) in enumerate(zip(LOSSES, LOSS_HUES, strict=True)):
            for index, (alpha, depth) in enumerate(zip(ALPHAS, ALPHA_DEPTHS, strict=True)):
                values = replicates[form, alpha, metric][field]
                score_bar(plot, position + (index - 1) * width, values, hue=hue, depth=depth, width=width)
                groups.append(values)
        setup_score_axis(plot, x, LOSSES.values(), groups, difference=difference, metric=metric)

    handles = []
    for alpha, depth in zip(ALPHAS, ALPHA_DEPTHS, strict=True):
        main, secondary = palette.color_pair(LOSS_HUES[0], depth)
        handles.append(Patch(facecolor=secondary, edgecolor=main, label=alpha))
    fig.legend(handles=handles, title="Delta coefficient α", ncols=3,
               loc="outside upper center", frameon=False)
    xlabel = "Delta loss metric"
    fig.supxlabel(xlabel if difference else f"{xlabel}\n{FOCUSED_Y_NOTE}")
    fig.supylabel("Paired score difference vs. MSE × 1" if difference else "Equal-fold score")
    plotting.save_figure(fig, output / FIGURES[field], dpi=180)
    plt.close(fig)


def render(data, output):
    scores = read_scores(data / "delta-metrics.csv")
    replicates = read_replicates(data / FOLD_SCORES, scores, coordinates=("loss", "alpha"),
                                 reference=("mse", "1"), folds=(0, 1, 2))
    for field in FIGURES:
        render_panels(replicates, field, output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == "__main__":
    main()
