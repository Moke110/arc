"""Render the Step 06.2 figures from the published score and epoch CSVs."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

STEP = Path(__file__).resolve().parents[2]
STAGE = STEP.parents[1]
sys.path[:0] = [str(STAGE / "util/output/api"), str(STEP / "execution")]

from architecture_roster import NO_X2, SHAPES, objective_label
from stage2_util import palette, plotting

LOWER_IS_BETTER = {"Δ RMSE", "Δ E-distance", "Δ MMD"}
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
SHAPE_COLORS = {
    "original": (20, 2),
    "output-delta": (60, 2),
    "output-delta-nox2": (60, 3),
    "z-interact_MLP": (140, 2),
    "z-interact_FiLM": (140, 3),
    "shared-z": (220, 2),
    "z-res": (260, 2),
    "shared-z-res": (300, 2),
}
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
    shapes = SHAPES
    folds = (0,) if smoke else range(3)
    if set(summary) != set(product(shapes, LABELS)) or set(replicates) != set(product(shapes, LABELS, folds)):
        raise ValueError("aggregate plotting CSV requires the complete architecture-by-metric grid")
    for (shape, metric), fields in summary.items():
        for field, mean in fields.items():
            values = [replicates[shape, metric, fold][field] for fold in folds]
            defined = [value for value in values if math.isfinite(value)]
            if defined:
                matches = math.isclose(mean, math.fsum(defined) / len(defined), rel_tol=0, abs_tol=1e-12)
            else:
                matches = math.isnan(mean)
            if not matches:
                raise ValueError(f"fold scores disagree with aggregate for {shape}, {metric}, {field}")
    return shapes, folds, summary, replicates


def render_metrics(data, output, smoke):
    shapes, folds, summary, replicates = read_scores(data, smoke)
    positions = range(len(shapes))
    for field, filename in FIELDS.items():
        difference = field == "score_difference"
        figure, plots = plt.subplots(3, 3, figsize=(18, 11), layout="constrained")
        for (metric, label), plot in zip(LABELS.items(), plots.flat, strict=True):
            plot.set_title(f"{label} {'↓' if metric in LOWER_IS_BETTER else '↑'}", fontweight="bold")
            means = [summary[shape, metric][field] for shape in shapes]
            groups = [[replicates[shape, metric, fold][field] for fold in folds] for shape in shapes]
            for position, shape, mean, values in zip(positions, shapes, means, groups, strict=True):
                hue, depth = SHAPE_COLORS[shape]
                plotting.replicate_bar(plot, position, values, hue=hue, depth=depth, width=.68)
                if math.isnan(mean):
                    plot.text(position, .02, "NA", ha="center", fontsize=7, transform=plot.get_xaxis_transform())
            reference = 0 if difference or metric == "Δ R²" else None
            lower, upper = plotting.focus_y(plot, groups, reference=reference)
            if lower <= 0 <= upper:
                plot.axhline(0, color="#666666", linewidth=.8, zorder=1.5)
            plot.set_xticks(positions, shapes, rotation=30, ha="right", fontsize=8)
            plot.set_xlim(-.5, len(shapes) - .5)
            plot.grid(axis="y", color="#777777", alpha=.25, zorder=1)
            plot.spines[["top", "right"]].set_visible(False)
        title = "Bars are means · Dots are replicates · Error bars are sample SD"
        if not difference:
            title += "\nY-axis ranges follow the data and may omit zero."
        figure.suptitle(title)
        figure.supylabel("Paired difference from the original architecture" if difference
                         else "Held-out score")
        plotting.save_figure(figure, output / filename, dpi=180)
        plt.close(figure)


def read_curves(path, smoke):
    curves = {shape: {} for shape in SHAPES}
    for row in plotting.read_rows(path):
        shape = row["shape"]
        if row["objective"] != objective_label(shape):
            raise ValueError(f"epoch CSV objective differs from the registered objective in {shape}")
        run = int(row["fold"]), int(row["seed"])
        points = curves[shape].setdefault(run, [])
        epoch = int(row["epoch"])
        values = tuple(float(row[field]) for field in ("train_loss", "validation_loss"))
        if epoch != len(points) + 1 or not all(math.isfinite(value) and value > 0 for value in values):
            raise ValueError(f"invalid objective in {shape}, run {run}, epoch {epoch}")
        points.append((epoch, *values))
    roster = set(next(iter(curves.values())))
    folds = (0,) if smoke else range(3)
    if roster != {(fold, 42) for fold in folds} or any(set(runs) != roster for runs in curves.values()):
        raise ValueError("epoch CSV requires the registered run roster across all architectures")
    return curves


def render_curves(data, output, smoke):
    curves = read_curves(data / "epoch-objectives.csv", smoke)
    figure, plots = plt.subplots(2, 4, figsize=(16, 7), sharex=True, layout="constrained")
    for (shape, runs), plot in zip(curves.items(), plots.flat, strict=True):
        if shape != NO_X2 and plot is not plots.flat[0]:
            plot.sharey(plots.flat[0])
        plotting.learning_curve(plot, [tuple(zip(*runs[run], strict=True)) for run in sorted(runs)])
        plot.set_title(f"{shape}\n{objective_label(shape)}", fontsize=10)
        plot.set_yscale("log")
        plot.grid(axis="y", color="#dedede", linewidth=.55)
        plot.spines[["top", "right"]].set_visible(False)
    handles = [Line2D([], [], color=palette.color_pair(hue).main, label=label)
               for hue, label in ((220, "Train"), (140, "Validation"))]
    figure.legend(handles=handles, ncols=len(handles), loc="outside upper center", frameon=False)
    figure.supxlabel("Epoch")
    figure.supylabel("Objective per condition (log scale)")
    plotting.save_figure(figure, output / CURVE_NAME, dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true", help="Render bounded fold-0 coordinates")
    args = parser.parse_args()
    render_metrics(args.data, args.output, args.smoke)
    render_curves(args.data, args.output, args.smoke)


if __name__ == "__main__":
    main()
