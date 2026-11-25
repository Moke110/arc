"""Render the Step 06.1 figures from the published score and epoch CSVs."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

STAGE = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(STAGE / "util/output/api"), str(STAGE / "05_tsg_loss/figure/script")]

from loss_plot_style import LOWER_IS_BETTER
from stage2_util import palette, plotting

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
MODIFICATIONS = {
    "Original": (220, ("widen_z-100",)),
    "Wider latent space": (140, ("widen_z-110", "widen_z-120", "widen_z-130")),
    "Wider encoder": (20, ("widen_enc",)),
    "Deeper network": (300, ("deepen_enc", "deepen_dec", "deepen-both")),
}
SHAPE_HUES = {shape: hue for hue, shapes in MODIFICATIONS.values() for shape in shapes}
FIELDS = {"value": "delta-metrics.png", "score_difference": "delta-metrics-paired-difference.png"}
CURVE_NAME = "learning-curves.png"


def score_table(path, replicate=False):
    summary = {}
    for row in plotting.read_rows(path):
        key = row["shape"], row["metric"]
        if replicate:
            key = (*key, int(row["fold"]))
        values = {field: float(row[field]) if row[field] else math.nan for field in FIELDS}
        if key in summary or any(math.isinf(value) for value in values.values()):
            raise ValueError(f"invalid or duplicate score in {path.name}: {key}")
        summary[key] = values
    return summary


def read_scores(data, smoke):
    summary = score_table(data / "delta-metrics.csv")
    replicates = score_table(data / "fold-scores.csv", replicate=True)
    shapes = tuple(shape for shape in SHAPE_HUES if not smoke or shape not in ("widen_z-110", "widen_z-120"))
    folds = (0,) if smoke else range(3)
    if set(summary) != set(product(shapes, LABELS)) or set(replicates) != set(product(shapes, LABELS, folds)):
        raise ValueError("aggregate plotting CSV requires the complete architecture-by-metric grid")
    return shapes, folds, summary, replicates


def render_metrics(data, output, smoke):
    shapes, folds, summary, replicates = read_scores(data, smoke)
    positions = range(len(shapes))
    for field, filename in FIELDS.items():
        difference = field == "score_difference"
        figure, plots = plt.subplots(3, 3, figsize=(15, 10), layout="constrained")
        for (metric, label), plot in zip(LABELS.items(), plots.flat, strict=True):
            plot.set_title(f"{label} {'↓' if metric in LOWER_IS_BETTER else '↑'}", fontweight="bold")
            means = [summary[shape, metric][field] for shape in shapes]
            if difference:
                plot.axhline(0, color="#666666", linewidth=.8, zorder=2)
            else:
                plotting.bar(plot, positions, means, hue=[SHAPE_HUES[shape] for shape in shapes],
                             width=.68, zorder=3)
            for position, shape, mean in zip(positions, shapes, means, strict=True):
                hue = SHAPE_HUES[shape]
                values = [replicates[shape, metric, fold][field] for fold in folds]
                center = position - .10 if difference else position
                plotting.replicate_dots(plot, center, values, hue=hue, s=26, zorder=4,
                                        edgecolors="white", linewidths=.55)
                if difference:
                    plotting.scatter(plot, [position + .25], [mean], hue=hue, marker="D", s=42,
                                     zorder=5, edgecolors="#222222", linewidths=.75)
                if math.isnan(mean):
                    plot.text(position, .02, "NA", ha="center", fontsize=7, transform=plot.get_xaxis_transform())
            plot.set_xticks(positions, shapes, rotation=30, ha="right", fontsize=8)
            plot.set_xlim(-.5, len(shapes) - .5)
            plot.grid(axis="y", color="#777777", alpha=.25, zorder=1)
            plot.spines[["top", "right"]].set_visible(False)
        handles = []
        for label, (hue, _) in MODIFICATIONS.items():
            color = palette.color_pair(hue, 2)
            handles.append(Patch(facecolor=color.secondary, edgecolor=color.main, label=label))
        handles.append(Line2D([], [], linestyle="none", marker="o", color="#666666",
                              markersize=5, label="Replicate"))
        handles.append(Line2D([], [], linestyle="none", marker="D", color="#666666",
                              markersize=5, label="Mean") if difference
                       else Patch(facecolor="#bbbbbb", edgecolor="#666666", label="Mean"))
        figure.legend(handles=handles, loc="outside upper center", ncols=len(handles), frameon=False)
        figure.supylabel("Paired difference from the original architecture" if difference
                         else "Held-out score")
        plotting.save_figure(figure, output / filename, dpi=180)
        plt.close(figure)


def read_curves(path):
    curves = {shape: [] for shape in SHAPE_HUES}
    for row in plotting.read_rows(path):
        shape = row["shape"]
        points = curves[shape]
        epoch = int(row["epoch"])
        values = tuple(float(row[field]) for field in ("mean_train_loss", "mean_validation_loss"))
        if epoch != len(points) + 1 or not all(math.isfinite(value) and value > 0 for value in values):
            raise ValueError(f"invalid aggregate objective in {shape} epoch {epoch}")
        points.append((epoch, *values))
    if any(not points for points in curves.values()):
        raise ValueError("epoch CSV requires a curve for every architecture")
    return curves


def render_curves(data, output):
    curves = read_curves(data / "epoch-objectives.csv")
    figure, plots = plt.subplots(3, 3, figsize=(13, 9), sharex=True, layout="constrained")
    for index, (shape, plot) in enumerate(zip(curves, plots.flat)):
        separate_scale = shape == "widen_enc"
        if index and not separate_scale:
            plot.sharey(plots.flat[0])
        if index % 3 and not separate_scale:
            plot.tick_params(labelleft=False)
        epochs, train, valid = zip(*curves[shape], strict=True)
        plotting.line(plot, epochs, train, hue=220, linestyle="--", linewidth=1.0)
        plotting.line(plot, epochs, valid, hue=220, linewidth=1.2)
        plot.set_title(f"{shape} (separate y-scale)" if separate_scale else shape, fontsize=10)
        plot.set_yscale("log")
        plot.grid(axis="y", color="#dedede", linewidth=.55)
        plot.spines[["top", "right"]].set_visible(False)
    plots.flat[-1].set_axis_off()
    plots.flat[5].tick_params(labelbottom=True)
    color = palette.color_pair(220).main
    handles = [Line2D([], [], color=color, linestyle=style, label=label)
               for label, style in (("Train", "--"), ("Validation", "-"))]
    figure.legend(handles=handles, ncols=len(handles), loc="outside upper center", frameon=False)
    figure.supxlabel("Epoch")
    figure.supylabel("Mean objective per condition (log scale)")
    plotting.save_figure(figure, output / CURVE_NAME, dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true", help="Render both metrics for the bounded smoke roster")
    args = parser.parse_args()
    render_metrics(args.data, args.output, args.smoke)
    if not args.smoke:
        render_curves(args.data, args.output)


if __name__ == "__main__":
    main()
