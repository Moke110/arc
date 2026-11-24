"""Render the TranSiGen axis comparison directly from the two report summary CSVs."""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP.parent / "util/output/api"))

from stage2_util.palette import color
from stage2_util.plotting import save_figure


def _axis_label(axis):
    family, size = axis.rsplit("-", 1)
    return f"{family.replace('hvg-deg', 'HVG-DEG').replace('denoise-deg', 'Denoise-DEG')}\n{size}"


def render_metrics(csv, output):
    data = pd.read_csv(csv)
    axes = data.axis.drop_duplicates().tolist()
    metrics = data.metric.drop_duplicates().tolist()
    if len(metrics) != 9 or data.duplicated(["axis", "metric"]).any() or len(data) != 9 * len(axes):
        raise ValueError("metric summary must have nine scores for every reported axis")
    fig, plots = plt.subplots(3, 3, figsize=(13.2, 10.8), layout="constrained")
    palette = [color("blue" if axis.startswith("hvg-") else "orange", 2) for axis in axes]
    for plot, metric in zip(plots.flat, metrics, strict=True):
        values = data.loc[data.metric.eq(metric)].set_index("axis").reindex(axes).value
        x = np.arange(len(axes))
        plot.bar(x, values, color=palette, width=.65)
        for index in np.flatnonzero(values.isna()):
            plot.text(index, .04, "NA", ha="center", transform=plot.get_xaxis_transform())
        plot.axhline(0, color="#777777", lw=.6)
        plot.set_title(metric)
        plot.set_xticks(x, [_axis_label(axis) for axis in axes], fontsize=8)
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
    image = plot.imshow(values, cmap="Blues", vmin=0, vmax=high, aspect="auto")
    plot.set_xticks(np.arange(len(axes)), [_axis_label(axis) for axis in axes], fontsize=9)
    plot.set_yticks(np.arange(len(diseases)), [name.replace("_", " ") for name in diseases])
    for (row, column), value in np.ndenumerate(values):
        plot.text(column, row, f"{value:.3f}" if np.isfinite(value) else "NA",
                  ha="center", va="center", color="white" if value > .66 * high else "#202020")
    plot.set_title("Signed disease DEG F1")
    fig.colorbar(image, ax=plot, label="F1")
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
