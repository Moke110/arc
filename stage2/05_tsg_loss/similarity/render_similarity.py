"""Render bounded Delta similarity sweep figures from the published CSVs."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LOSS = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(LOSS.parent / "util/output/api"), str(LOSS / "figure/script"),
                str(LOSS / "05.3_mse-l1_mse-x/figure/script")]

from loss_plot_style import FOLD_HUES, curve_legend_handles, metric_direction
from render_metrics import FIGURES, LABELS
from render_learning_curves import EPOCH_SOURCE, FIGURE_NAME, FOLDS, read_curves
from stage2_util import plotting

# These are the accepted CSV report coordinates. Importing them from the
# objective would require PyTorch for CSV-only rendering. Moving them out
# of the objective would change the code digest sealed by run registrations.
ALPHAS = ("0", "0.1", "0.2", "0.3", "0.4", "0.5")
REFERENCE_ALPHA = "0"
SMOKE_ALPHAS = (REFERENCE_ALPHA, "0.5")
METRIC_LABELS = {"pcc": "PCC", "cos": "cosine"}



def read_scores(path, alphas):
    scores = {}
    for row in plotting.read_rows(path):
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
        hue = metric_direction(plot, metric, label)
        plotting.line(plot, x, values, hue=hue, marker="o", markersize=3.5, zorder=3)
        for position, value in zip(x, values, strict=True):
            if math.isnan(value):
                plot.text(position, .025, "NA", ha="center", fontsize=7,
                          transform=plot.get_xaxis_transform(), zorder=4)
        if difference:
            limit = max(map(abs, plot.get_ylim()))
            plot.set_ylim(-limit, limit)
        plot.set(ylim=plot.get_ylim(), xticks=x, xticklabels=alphas)
        if len(alphas) > len(SMOKE_ALPHAS):
            for tick in plot.get_xticklabels():
                tick.set(rotation=45, ha="right", fontsize=9)
        plot.axhline(0, color="#666666", linewidth=.8, zorder=2)
        plot.grid(axis="y", color="#777777", alpha=.25, zorder=1)
        plot.spines[["top", "right"]].set_visible(False)

    fig.supxlabel(f"Delta {METRIC_LABELS[similarity]} loss weight α")
    fig.supylabel("Signed score difference vs. α = 0" if difference else "Equal-fold score")
    plotting.save_figure(fig, output / FIGURES[field], dpi=180)
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
            plotting.line(plot, epochs, train, hue=FOLD_HUES[fold], linestyle="--", linewidth=1.0)
            plotting.line(plot, epochs, valid, hue=FOLD_HUES[fold], linewidth=1.2)
            best = curves[alpha, fold][0][0][2]
            plotting.scatter(plot, best, valid[best - 1], hue=FOLD_HUES[fold], s=12, zorder=3)
        plot.set_title(f"α = {alpha}", fontsize=10)
        plot.set_yscale("log")
        plot.grid(axis="y", color="#dedede", linewidth=.55)
        plot.spines[["top", "right"]].set_visible(False)

    handles = curve_legend_handles(folds, selected_epoch=True)
    fig.legend(handles=handles, ncols=len(handles), loc="outside upper center", frameon=False)
    fig.supxlabel(f"Epoch · Delta {METRIC_LABELS[metric]} loss")
    fig.supylabel("Objective per condition (log scale)")
    plotting.save_figure(fig, output / FIGURE_NAME, dpi=180)
    plt.close(fig)



def render(data, output, metric, *, smoke=False):
    if metric not in METRIC_LABELS:
        raise ValueError(f"unsupported Delta similarity metric: {metric}")
    alphas = SMOKE_ALPHAS if smoke else ALPHAS
    scores = read_scores(data / "delta-metrics.csv", alphas)
    for field in FIGURES:
        render_panels(scores, field, output, alphas, metric)
    render_curves(data, output, metric, smoke=smoke)


def main(metric=None):
    parser = argparse.ArgumentParser(description=__doc__)
    if metric is None:
        parser.add_argument("--metric", choices=METRIC_LABELS, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    render(args.data, args.output, metric or args.metric, smoke=args.smoke)


if __name__ == "__main__":
    main()
