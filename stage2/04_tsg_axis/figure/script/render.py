"""Render the TranSiGen axis comparison from published summary and replicate CSVs."""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FIGURE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FIGURE.parents[1] / "util/output/api"))

from stage2_util.plotting import replicate_bar, save_figure

LOWER_IS_BETTER = {"Δ RMSE", "Δ E-distance", "Δ MMD"}
AXIS_HUES = {"hvg": 220, "denoise": 20}


def panels(summary_csv, fold_csv, category):
    summary = pd.read_csv(summary_csv)
    folds = pd.read_csv(fold_csv)
    axes = summary.axis.drop_duplicates().tolist()
    categories = summary[category].drop_duplicates().tolist()
    groups = ["axis", category]
    if summary.duplicated(groups).any() or len(summary) != len(axes) * len(categories):
        raise ValueError(f"summary must contain each axis and {category}")
    replicate_keys = [*groups, "fold"]
    if "seed" in folds:
        replicate_keys.append("seed")
    if folds.duplicated(replicate_keys).any():
        raise ValueError(f"duplicate replicate in {fold_csv}")

    replicates = folds.groupby(groups)["value"]
    means = replicates.agg(["size", "count", "mean"]).reindex(summary.set_index(groups).index)
    if (len(replicates) != len(summary)
            or not np.array_equal(means["size"], summary.n_folds)
            or not np.array_equal(means["count"], summary.n_eligible_folds)
            or not np.allclose(means["mean"], summary.value, rtol=0, atol=1e-12, equal_nan=True)):
        raise ValueError(f"replicates in {fold_csv} disagree with {summary_csv}")
    return axes, categories, replicates


def render_panel(plot, axes, replicates, category):
    for x, axis in enumerate(axes):
        method, *_, size = axis.split("-")
        values = replicates.get_group((axis, category))
        if np.isfinite(values).any():
            replicate_bar(plot, x, values, hue=AXIS_HUES[method],
                          depth=5 - int(size.removesuffix("k")), width=.66)
        else:
            plot.text(x, .04, "NA", ha="center", transform=plot.get_xaxis_transform())
    plot.axhline(0, color="#777777", lw=.6)
    plot.set_xticks(np.arange(len(axes)), axes, fontsize=8, rotation=30,
                    ha="right", rotation_mode="anchor")
    plot.spines[["top", "right"]].set_visible(False)
    plot.grid(axis="y", color="#eeeeee")
    plot.set_axisbelow(True)


def render_metrics(csv, fold_csv, output):
    axes, metrics, replicates = panels(csv, fold_csv, "metric")
    if len(metrics) != 9:
        raise ValueError("metric summary must have nine scores for every reported axis")
    fig, plots = plt.subplots(3, 3, figsize=(13.2, 10.8), layout="constrained")
    for plot, metric in zip(plots.flat, metrics, strict=True):
        render_panel(plot, axes, replicates, metric)
        direction = "↓" if metric in LOWER_IS_BETTER else "↑"
        plot.set_title(f"{metric} {direction}")
    fig.suptitle("TranSiGen Delta scores on fixed tenv axes · mean ± SD")
    save_figure(fig, output / "delta-metrics.png", dpi=180)
    plt.close(fig)


def render_disease(csv, fold_csv, output):
    axes, diseases, replicates = panels(csv, fold_csv, "disease")
    fig, plots = plt.subplots((len(diseases) + 1) // 2, 2, figsize=(11, 7),
                              sharey=True, squeeze=False, layout="constrained")
    for plot, disease in zip(plots.flat, diseases):
        render_panel(plot, axes, replicates, disease)
        plot.set_title(disease.replace("_", " "))
        plot.set_ylabel("F1")
    for plot in list(plots.flat)[len(diseases):]:
        plot.remove()
    fig.suptitle("Signed disease DEG F1 · mean ± SD")
    save_figure(fig, output / "disease-f1.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=FIGURE / "data")
    parser.add_argument("--fold-data", type=Path, help="directory with replicate CSVs, if separate")
    parser.add_argument("--output", type=Path, default=FIGURE)
    args = parser.parse_args()
    fold_data = args.fold_data or args.data
    render_metrics(args.data / "delta-metrics.csv", fold_data / "fold-metrics.csv", args.output)
    render_disease(args.data / "disease-f1.csv", fold_data / "fold-disease-f1.csv", args.output)


if __name__ == "__main__":
    main()
