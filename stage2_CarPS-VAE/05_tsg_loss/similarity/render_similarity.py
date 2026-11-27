"""Render bounded Delta similarity sweep figures from the published CSVs."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOSS = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(LOSS.parent / "util/output/api"), str(LOSS / "execution"),
                str(LOSS / "05.3_mse-l1_mse-x/execution")]

from loss_plot_style import FOLD_SCORES, read_replicates, render_weight_panels
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

def render_curves(data, output, metric, *, smoke):
    alphas = SMOKE_ALPHAS if smoke else ALPHAS
    folds = (0,) if smoke else FOLDS
    curves = read_curves(data / EPOCH_SOURCE, alphas, folds, 2 if smoke else 1000)
    shape = (1, 2) if smoke else (2, 3)
    fig, plots = plt.subplots(*shape, figsize=(9, 4) if smoke else (12, 8),
                              sharex=True, sharey=True, layout="constrained")
    for alpha, plot in zip(alphas, plots.flat, strict=True):
        runs = []
        selected_epochs = []
        for fold in folds:
            points = curves[alpha, fold]
            _, epochs, train, valid, _ = zip(*points, strict=True)
            runs.append((epochs, train, valid))
            best = points[0][0][2]
            selected_epochs.append((best, valid[best - 1]))
        plotting.learning_curve(plot, runs)
        best_epochs, best_losses = zip(*selected_epochs, strict=True)
        plotting.scatter(plot, best_epochs, best_losses, hue=140, s=12, zorder=3,
                         label="Selected epoch")
        plot.set_title(f"α = {alpha}", fontsize=10)
        plot.set_yscale("log")
        plot.grid(axis="y", color="#dedede", linewidth=.55)
        plot.spines[["top", "right"]].set_visible(False)

    handles, labels = plots.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncols=3, loc="outside upper center", frameon=False)
    fig.supxlabel(f"Epoch · Delta {METRIC_LABELS[metric]} loss")
    fig.supylabel("Objective per condition (log scale)")
    plotting.save_figure(fig, output / FIGURE_NAME, dpi=180)
    plt.close(fig)


def render(data, output, metric, *, smoke=False):
    if metric not in METRIC_LABELS:
        raise ValueError(f"unsupported Delta similarity metric: {metric}")
    alphas = SMOKE_ALPHAS if smoke else ALPHAS
    scores = read_scores(data / "delta-metrics.csv", alphas)
    replicates = read_replicates(data / FOLD_SCORES, scores, coordinates=("alpha",),
                                 reference=(REFERENCE_ALPHA,), folds=(0,) if smoke else FOLDS)
    for field, filename in FIGURES.items():
        render_weight_panels(
            replicates, field, output / filename, alphas, LABELS,
            xlabel=f"Delta {METRIC_LABELS[metric]} loss weight α",
            reference=REFERENCE_ALPHA, rotate=not smoke,
        )
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
