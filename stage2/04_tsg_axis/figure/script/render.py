"""Render the TranSiGen axis comparison directly from the two report summary CSVs."""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.image import BboxImage
import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP.parent / "util/output/api"))

from stage2_util.palette import color_pair
from stage2_util.plotting import bar, heatmap, save_figure

LOWER_IS_BETTER = {"Δ RMSE", "Δ E-distance", "Δ MMD"}
AXIS_HUES = {"hvg": 220, "denoise": 20}


def _add_metric_background(plot, metric):
    hue = 340 if metric in LOWER_IS_BETTER else 220
    gradient = np.empty((256, 1, 4))
    gradient[:, :, :3] = color_pair(hue).main[:3]
    gradient[:, 0, 3] = np.linspace(0, .25, len(gradient))
    background = BboxImage(plot.bbox, origin="lower", interpolation="bilinear", zorder=0)
    background.set_data(gradient)
    plot.add_artist(background)


def render_metrics(csv, output):
    data = pd.read_csv(csv)
    axes = data.axis.drop_duplicates().tolist()
    metrics = data.metric.drop_duplicates().tolist()
    if len(metrics) != 9 or data.duplicated(["axis", "metric"]).any() or len(data) != 9 * len(axes):
        raise ValueError("metric summary must have nine scores for every reported axis")
    fig, plots = plt.subplots(3, 3, figsize=(13.2, 10.8), layout="constrained")
    hues = [AXIS_HUES[axis.split("-")[0]] for axis in axes]
    depths = [5 - int(axis.rsplit("-", 1)[-1].removesuffix("k")) for axis in axes]
    for plot, metric in zip(plots.flat, metrics, strict=True):
        values = data.loc[data.metric.eq(metric)].set_index("axis").reindex(axes).value
        x = np.arange(len(axes))
        _add_metric_background(plot, metric)
        bar(plot, x, values, hue=hues, depth=depths, width=.65)
        for index in np.flatnonzero(values.isna()):
            plot.text(index, .04, "NA", ha="center", transform=plot.get_xaxis_transform())
        plot.axhline(0, color="#777777", lw=.6)
        direction = "↓" if metric in LOWER_IS_BETTER else "↑"
        plot.set_title(f"{metric} {direction}")
        plot.set_xticks(x, axes, fontsize=8, rotation=30, ha="right", rotation_mode="anchor")
        plot.spines[["top", "right"]].set_visible(False)
    fig.suptitle("TranSiGen Delta scores on fixed tenv axes")
    save_figure(fig, output / "delta-metrics.png", dpi=180)
    plt.close(fig)


def render_disease(csv, output):
    data = pd.read_csv(csv)
    axes = data.axis.drop_duplicates().tolist()
    diseases = data.disease.drop_duplicates().tolist()
    if data.duplicated(["axis", "disease"]).any() or len(data) != len(axes) * len(diseases):
        raise ValueError("disease summary must contain each axis and disease")
    values = data.pivot(index="disease", columns="axis", values="value").loc[diseases, axes].to_numpy()
    finite = values[np.isfinite(values)]
    high = float(np.max(finite)) if finite.size and np.max(finite) > 0 else 1.0
    fig, plot = plt.subplots(figsize=(max(5.5, 1.2 * len(axes) + 2.3), 4.4), layout="constrained")
    heatmap(plot, values, hue=220, vmin=0, vmax=high, aspect="auto", colorbar_label="F1")
    plot.set_xticks(np.arange(len(axes)), axes, fontsize=9)
    plot.set_yticks(np.arange(len(diseases)), [name.replace("_", " ") for name in diseases])
    for (row, column), value in np.ndenumerate(values):
        plot.text(column, row, f"{value:.3f}" if np.isfinite(value) else "NA",
                  ha="center", va="center", color="white" if value >= .40 * high else "#111111")
    plot.set_title("Signed disease DEG F1")
    save_figure(fig, output / "disease-f1.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render_metrics(args.data / "delta-metrics.csv", args.output)
    render_disease(args.data / "disease-f1.csv", args.output)


if __name__ == "__main__":
    main()
