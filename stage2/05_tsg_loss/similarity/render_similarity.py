"""Render bounded Delta similarity sweep figures from the published CSVs."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np

from similarity_objective import ALPHAS, METRICS, METRIC_LABELS, REFERENCE_ALPHA, SMOKE_ALPHAS

RECONSTRUCTION_FIGURES = Path(__file__).resolve().parents[1] / "05.3_mse-l1_mse-x/figure/script"
sys.path.append(str(RECONSTRUCTION_FIGURES))

from render_metrics import FIGURES, LABELS, LOWER_IS_BETTER
from render_learning_curves import COLORS, EPOCH_SOURCE, FIGURE_NAME, FOLDS, read_curves
from stage2_util.palette import color
from stage2_util.plotting import read_rows, save_figure


def read_scores(path, alphas):
    scores = {}
    for row in read_rows(path):
        key = row["alpha"], row["metric"]
        if key in scores:
            raise ValueError(f"duplicate score in complete sweep: {key}")
        record = {name: float(row[name]) if row[name] else math.nan
                  for name in (*FIGURES, "reference_weight_score")}
        if any(math.isinf(value) for value in record.values()):
            raise ValueError(f"infinite score for {key}")
        scores[key] = record
    if set(scores) != set(product(alphas, LABELS)):
        raise ValueError("metric summary must contain the complete sweep")
    for key, record in scores.items():
        baseline = scores[REFERENCE_ALPHA, key[1]]["value"]
        actual = record["reference_weight_score"], record["score_difference"]
        expected = baseline, record["value"] - baseline
        if any(a != b and not (math.isnan(a) and math.isnan(b))
               for a, b in zip(actual, expected, strict=True)):
            raise ValueError(f"signed score difference disagrees with the zero-weight score: {key}")
    return scores


def render_panels(scores, field, output, alphas, similarity):
    difference = field == "score_difference"
    fig, plots = plt.subplots(3, 3, figsize=(13, 9), layout="constrained")
    x = np.asarray(alphas, dtype=float)
    for (metric, label), plot in zip(LABELS.items(), plots.flat, strict=True):
        values = [scores[alpha, metric][field] for alpha in alphas]
        hue = "magenta" if metric in LOWER_IS_BETTER else "blue"
        plot.plot(x, values, color=color(hue, 2), marker="o", markersize=3.5, zorder=3)
        for position, value in zip(x, values, strict=True):
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
        plot.set(xlim=limits[0], ylim=limits[1], xticks=x, xticklabels=alphas, title=label)
        if len(alphas) > len(SMOKE_ALPHAS):
            for tick in plot.get_xticklabels():
                tick.set(rotation=45, ha="right", fontsize=9)
        plot.title.set_fontweight("bold")
        plot.axhline(0, color="#666666", linewidth=.8, zorder=2)
        plot.grid(axis="y", color="#777777", alpha=.25, zorder=1)
        plot.spines[["top", "right"]].set_visible(False)

    fig.supxlabel(f"Delta {METRIC_LABELS[similarity]} loss weight α")
    fig.supylabel("Signed score difference vs. α = 0" if difference else "Equal-fold score")
    save_figure(fig, output / FIGURES[field], dpi=180)
    plt.close(fig)


def render_curves(data, output, metric, *, smoke):
    alphas = SMOKE_ALPHAS if smoke else ALPHAS
    folds = (0,) if smoke else FOLDS
    curves = read_curves(data / EPOCH_SOURCE, alphas, folds, 2 if smoke else 1000)
    shape = (1, 2) if smoke else (2, 3)
    fig, plots = plt.subplots(*shape, figsize=(9, 4) if smoke else (12, 8),
                              sharex=True, sharey=True, layout="constrained")
    for alpha, plot in zip(alphas, plots.flat, strict=True):
        for fold in folds:
            _, epochs, train, valid, _ = zip(*curves[alpha, fold], strict=True)
            plot.plot(epochs, train, linestyle="--", color=COLORS[fold], linewidth=1.0)
            plot.plot(epochs, valid, color=COLORS[fold], linewidth=1.2)
            best = curves[alpha, fold][0][0][2]
            plot.scatter(best, valid[best - 1], color=COLORS[fold], s=12, zorder=3)
        plot.set_title(f"α = {alpha}", fontsize=10)
        plot.set_yscale("log")
        plot.grid(axis="y", color="#dedede", linewidth=.55)
        plot.spines[["top", "right"]].set_visible(False)

    handles = [Line2D([], [], color=COLORS[fold], label=f"Fold {fold}") for fold in folds]
    handles.extend((Line2D([], [], color="#333333", linestyle="--", label="Train"),
                    Line2D([], [], color="#333333", label="Validation"),
                    Line2D([], [], color="#333333", linestyle="None", marker="o", label="Selected epoch")))
    fig.legend(handles=handles, ncols=len(handles), loc="outside upper center", frameon=False)
    fig.supxlabel(f"Epoch · Delta {METRIC_LABELS[metric]} loss")
    fig.supylabel("Objective per condition (log scale)")
    save_figure(fig, output / FIGURE_NAME, dpi=180)
    plt.close(fig)



def render(data, output, metric, *, smoke=False):
    if metric not in METRICS:
        raise ValueError(f"unsupported Delta similarity metric: {metric}")
    alphas = SMOKE_ALPHAS if smoke else ALPHAS
    scores = read_scores(data / "delta-metrics.csv", alphas)
    for field in FIGURES:
        render_panels(scores, field, output, alphas, metric)
    render_curves(data, output, metric, smoke=smoke)


def main(metric=None):
    parser = argparse.ArgumentParser(description=__doc__)
    if metric is None:
        parser.add_argument("--metric", choices=METRICS, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    render(args.data, args.output, metric or args.metric, smoke=args.smoke)


if __name__ == "__main__":
    main()
