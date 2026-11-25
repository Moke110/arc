"""Render equal-fold Delta scores and signed differences against native MSE."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LOSS = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(LOSS.parent / "util/output/api"), str(LOSS / "figure/script")]

from loss_plot_style import metric_direction
from stage2_util import plotting

LOSSES = {"mse": "MSE", "cos": "Cos", "pcc": "PCC", "l1": "L1"}
ALPHAS = ("0.75", "1", "1.25")
ALPHA_HUES = (220, 20, 140)
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


def render_panels(scores, field, output):
    difference = field == "score_difference"
    fig, plots = plt.subplots(3, 3, figsize=(13, 9), layout="constrained")
    x = np.arange(len(LOSSES))
    for (metric, label), plot in zip(LABELS.items(), plots.flat, strict=True):
        metric_direction(plot, metric, label)
        groups = [([scores[form, alpha, metric][field] for form in LOSSES], hue, 2, alpha)
                  for alpha, hue in zip(ALPHAS, ALPHA_HUES, strict=True)]
        bars = plotting.grouped_bar(plot, x, groups, width=.72, linewidth=.55, zorder=3)
        for group_bars in bars:
            for bar in group_bars:
                if math.isnan(bar.get_height()):
                    plot.text(bar.get_center()[0], .025, "NA", ha="center", fontsize=7,
                              transform=plot.get_xaxis_transform(), zorder=4)
        if difference:
            limit = max(map(abs, plot.get_ylim()))
            plot.set_ylim(-limit, limit)
        plot.set(ylim=plot.get_ylim(), xticks=x, xticklabels=LOSSES.values())
        plot.axhline(0, color="#666666", linewidth=.8, zorder=2)
        plot.grid(axis="y", color="#777777", alpha=.25, zorder=1)
        plot.spines[["top", "right"]].set_visible(False)

    handles, labels = plots[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels,
               title="Delta coefficient α   ·   ↑ Higher is better   ·   ↓ Lower is better",
               ncols=3, loc="outside upper center", frameon=False)
    fig.supxlabel("Delta loss metric")
    fig.supylabel("Score difference vs. MSE × 1" if difference else "Equal-fold score")
    plotting.save_figure(fig, output / FIGURES[field], dpi=180)
    plt.close(fig)


def render(data, output):
    scores = read_scores(data / "delta-metrics.csv")
    for field in FIGURES:
        render_panels(scores, field, output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == "__main__":
    main()
