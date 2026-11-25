"""Render the four Delta-metric panels from the equal-fold scores CSV."""

import argparse
import csv
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP.parent / "util/output/api"))

from stage2_util.palette import color, entry_hue
from stage2_util.plotting import save_figure

LOSSES = {"mse": "MSE", "cos": "Cos", "pcc": "PCC", "l1": "L1"}
ALPHAS = ("1/4", "1/3", "1/2")
LABELS = {
    "JAC_top200": "JAC\ntop200",
    "JAC_pos100": "JAC\npos100",
    "JAC_neg100": "JAC\nneg100",
    "Δ PCC": "Δ PCC",
    "Δ Spearman": "Δ Spearman",
    "Δ RMSE": "Δ RMSE",
    "Δ E-distance": "Δ E-distance",
    "Δ MMD": "Δ MMD",
    "Δ R²": "Δ R²",
}


def read_scores(path):
    scores = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            key = row["loss"], row["alpha"], row["metric"]
            if key in scores:
                raise ValueError(f"duplicate score in complete loss grid: {key}")
            value = float(row["value"]) if row["value"] else math.nan
            if math.isinf(value):
                raise ValueError(f"infinite score for {key}")
            scores[key] = value
    if set(scores) != set(product(LOSSES, ALPHAS, LABELS)):
        raise ValueError("metric summary must contain the complete loss grid")
    return scores


def render(data, output):
    scores = read_scores(data / "delta-metrics.csv")
    fig, plots = plt.subplots(2, 2, figsize=(14, 8), sharey=True, layout="constrained")
    x = np.arange(len(LABELS))
    for (form, label), plot in zip(LOSSES.items(), plots.flat, strict=True):
        for index, alpha in enumerate(ALPHAS):
            positions = x + (index - 1) * .25
            values = [scores[form, alpha, metric] for metric in LABELS]
            plot.bar(positions, values, color=color(entry_hue(index), 5), width=.23,
                     edgecolor=color(entry_hue(index), 2), linewidth=.4, label=alpha)
            for position, value in zip(positions, values, strict=True):
                if math.isnan(value):
                    plot.text(position, .025, "NA", ha="center", fontsize=7,
                              transform=plot.get_xaxis_transform())
        plot.set_xticks(x, LABELS.values(), fontsize=8)
        plot.set_title(label)
        plot.axhline(0, color="#777777", linewidth=.65)
        plot.spines[["top", "right"]].set_visible(False)
    handles, labels = plots[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, title="Delta loss weight α", ncols=3,
               loc="outside upper center", frameon=False)
    fig.supylabel("Equal-fold score")
    save_figure(fig, output / "delta-metrics.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == "__main__":
    main()
