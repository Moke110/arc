"""Render equal-fold Delta scores and signed differences against native MSE."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "stage2/util/output/api"))

from stage2_util.palette import color
from stage2_util.plotting import read_rows, save_figure

LOSSES = {"mse": "MSE", "cos": "Cos", "pcc": "PCC", "l1": "L1"}
ALPHAS = ("0.75", "1", "1.25")
COLORS = ("blue", "orange", "green")
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
LOWER_IS_BETTER = {"Δ RMSE", "Δ E-distance", "Δ MMD"}
FIGURES = {"value": "delta-metrics.png",
           "score_difference": "delta-metrics-improvement.png"}


def read_scores(path):
    scores = {}
    for row in read_rows(path):
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
        for index, (alpha, hue) in enumerate(zip(ALPHAS, COLORS, strict=True)):
            positions = x + (index - 1) * .24
            values = [scores[form, alpha, metric][field] for form in LOSSES]
            plot.bar(positions, values, color=color(hue, 5), width=.22,
                     edgecolor=color(hue, 2), linewidth=.45, label=alpha, zorder=3)
            for position, value in zip(positions, values, strict=True):
                if math.isnan(value):
                    plot.text(position, .025, "NA", ha="center", fontsize=7,
                              transform=plot.get_xaxis_transform(), zorder=4)
        if difference:
            limit = max(map(abs, plot.get_ylim()))
            plot.set_ylim(-limit, limit)
        limits = plot.get_xlim(), plot.get_ylim()
        tint = "#efc5c5" if metric in LOWER_IS_BETTER else "#c1d0e5"
        gradient = LinearSegmentedColormap.from_list("direction", ["#fdfdfd", tint])
        plot.imshow([[0.], [1.]], cmap=gradient, origin="lower", aspect="auto",
                    extent=(0, 1, 0, 1), transform=plot.transAxes,
                    interpolation="bicubic", zorder=0)
        plot.set(xlim=limits[0], ylim=limits[1], xticks=x, xticklabels=LOSSES.values(), title=label)
        plot.title.set_fontweight("bold")
        plot.axhline(0, color="#666666", linewidth=.8, zorder=2)
        plot.grid(axis="y", color="#777777", alpha=.25, zorder=1)
        plot.spines[["top", "right"]].set_visible(False)

    handles, labels = plots[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels,
               title="Delta coefficient α   ·   Blue panels favor higher scores   ·   Red panels favor lower scores",
               ncols=3, loc="outside upper center", frameon=False)
    fig.supxlabel("Delta loss metric")
    fig.supylabel("Score difference vs. MSE × 1" if difference else "Equal-fold score")
    save_figure(fig, output / FIGURES[field], dpi=180)
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
