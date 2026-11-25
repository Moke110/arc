"""Render equal-fold Delta scores and signed differences against zero Delta L1 weight."""

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

ALPHAS = ("0", "0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75", "2")
SMOKE_ALPHAS = ("0", "1")
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


def read_scores(path, alphas):
    scores = {}
    for row in plotting.read_rows(path):
        key = row["alpha"], row["metric"]
        if key in scores:
            raise ValueError(f"duplicate score in complete sweep: {key}")
        record = {name: float(row[name]) if row[name] else math.nan
                  for name in (*FIGURES, "zero_weight_score")}
        if any(math.isinf(value) for value in record.values()):
            raise ValueError(f"infinite score for {key}")
        scores[key] = record
    if set(scores) != set(product(alphas, LABELS)):
        raise ValueError("metric summary must contain the complete sweep")
    for key, record in scores.items():
        baseline = scores["0", key[1]]["value"]
        actual = record["zero_weight_score"], record["score_difference"]
        expected = baseline, record["value"] - baseline
        if any(a != b and not (math.isnan(a) and math.isnan(b))
               for a, b in zip(actual, expected, strict=True)):
            raise ValueError(f"signed score difference disagrees with the zero-weight score: {key}")
    return scores


def render_panels(scores, field, output, alphas):
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
        plot.axhline(0, color="#666666", linewidth=.8, zorder=2)
        plot.grid(axis="y", color="#777777", alpha=.25, zorder=1)
        plot.spines[["top", "right"]].set_visible(False)

    fig.supxlabel("Additive Delta L1 weight α")
    fig.supylabel("Signed score difference vs. α = 0" if difference else "Equal-fold score")
    plotting.save_figure(fig, output / FIGURES[field], dpi=180)
    plt.close(fig)


def render(data, output, *, smoke=False):
    alphas = SMOKE_ALPHAS if smoke else ALPHAS
    scores = read_scores(data / "delta-metrics.csv", alphas)
    for field in FIGURES:
        render_panels(scores, field, output, alphas)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    render(args.data, args.output, smoke=args.smoke)


if __name__ == "__main__":
    main()
