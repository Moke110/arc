"""Render the Step 06.3 figures from the published score and epoch CSVs."""

import argparse
from itertools import product
import math
import statistics
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

STEP = Path(__file__).resolve().parents[2]
STAGE = STEP.parents[1]
sys.path[:0] = [str(STAGE / "util/output/api"), str(STEP / "execution")]

from composition_roster import SHAPES, coordinates, objective_label
from composition_scores import CONTRASTS, PRIMARY_METRICS
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
SHAPE_COLORS = dict(zip(SHAPES, ((20, 2), (60, 2), (60, 3), (140, 2), (220, 2), (220, 3)), strict=True))
FIELDS = {"value": "delta-metrics.png", "score_difference": "delta-metrics-paired-difference.png"}
CURVE_NAME = "learning-curves.png"
CONTRAST_NAME = "primary-contrasts.png"


def replicate_coordinates(smoke):
    return tuple(dict.fromkeys((fold, seed) for _, fold, seed in coordinates(smoke)))


def number(value):
    return float(value) if value else math.nan


def mean_matches(mean, values):
    defined = [value for value in values if math.isfinite(value)]
    return math.isclose(mean, statistics.mean(defined), rel_tol=0, abs_tol=1e-12) if defined else math.isnan(mean)


def score_table(path, replicate=False):
    summary = {}
    for row in plotting.read_rows(path):
        key = row["shape"], row["metric"]
        if replicate:
            key = (*key, int(row["fold"]), int(row["seed"]))
        values = {field: number(row[field]) for field in FIELDS}
        if key in summary or any(math.isinf(value) for value in values.values()):
            raise ValueError(f"invalid or duplicate score in {path.name}: {key}")
        summary[key] = values
    return summary


def read_scores(data, smoke):
    summary = score_table(data / "delta-metrics.csv")
    replicates = score_table(data / "fold-scores.csv", replicate=True)
    keys = replicate_coordinates(smoke)
    if (set(summary) != set(product(SHAPES, LABELS))
            or set(replicates) != {(shape, metric, *key) for shape, metric, key in product(SHAPES, LABELS, keys)}):
        raise ValueError("plotting CSV requires all registered arm, metric, fold and seed rows")
    for (shape, metric), fields in summary.items():
        for field, mean in fields.items():
            values = [replicates[shape, metric, *key][field] for key in keys]
            if not mean_matches(mean, values):
                raise ValueError(f"coordinate scores disagree with aggregate for {shape}, {metric}, {field}")
    return keys, summary, replicates


def render_metrics(data, output, smoke):
    keys, summary, replicates = read_scores(data, smoke)
    positions = range(len(SHAPES))
    for field, filename in FIELDS.items():
        difference = field == "score_difference"
        figure, plots = plt.subplots(3, 3, figsize=(18, 11), layout="constrained")
        for (metric, label), plot in zip(LABELS.items(), plots.flat, strict=True):
            plot.set_title(f"{label} {'↓' if metric in LOWER_IS_BETTER else '↑'}", fontweight="bold")
            means = [summary[shape, metric][field] for shape in SHAPES]
            groups = [[replicates[shape, metric, *key][field] for key in keys] for shape in SHAPES]
            for position, shape, mean, values in zip(positions, SHAPES, means, groups, strict=True):
                hue, depth = SHAPE_COLORS[shape]
                plotting.replicate_bar(plot, position, values, hue=hue, depth=depth, width=.68)
                if math.isnan(mean):
                    plot.text(position, .02, "NA", ha="center", fontsize=7, transform=plot.get_xaxis_transform())
            reference = 0 if difference or metric == "Δ R²" else None
            lower, upper = plotting.focus_y(plot, groups, reference=reference)
            if lower <= 0 <= upper:
                plot.axhline(0, color="#666666", linewidth=.8, zorder=1.5)
            plot.set_xticks(positions, SHAPES, rotation=30, ha="right", fontsize=8)
            plot.set_xlim(-.5, len(SHAPES) - .5)
            plot.grid(axis="y", color="#777777", alpha=.25, zorder=1)
            plot.spines[["top", "right"]].set_visible(False)
        title = "Bars are means · Dots are fold and seed coordinates · Error bars are sample SD"
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
    roster = set(replicate_coordinates(smoke))
    if any(set(runs) != roster for runs in curves.values()):
        raise ValueError("epoch CSV requires the registered run roster across all architectures")
    return curves


def render_curves(data, output, smoke):
    curves = read_curves(data / "epoch-objectives.csv", smoke)
    figure, plots = plt.subplots(2, 3, figsize=(13, 8), sharex=True, layout="constrained")
    family_axes = {}
    for (shape, runs), plot in zip(curves.items(), plots.flat, strict=True):
        reference = family_axes.setdefault(objective_label(shape), plot)
        if plot is not reference:
            plot.sharey(reference)
        plotting.learning_curve(plot, [tuple(zip(*runs[run], strict=True)) for run in replicate_coordinates(smoke)])
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


def contrast_values(path, replicate):
    values = {}
    for row in plotting.read_rows(path):
        key = row["contrast"], row["metric"]
        if replicate:
            key = (*key, int(row["fold"]), int(row["seed"]))
        value = number(row["value"])
        if key in values or math.isinf(value):
            raise ValueError(f"invalid or duplicate contrast in {path.name}: {key}")
        values[key] = value
    return values


def render_contrasts(data, output, smoke):
    keys = replicate_coordinates(smoke)
    replicates = contrast_values(data / "primary-contrasts.csv", replicate=True)
    summary = contrast_values(data / "primary-contrast-summary.csv", replicate=False)
    if (set(summary) != set(product(CONTRASTS, PRIMARY_METRICS))
            or set(replicates) != {(contrast, metric, *key)
                                  for contrast, metric, key in product(CONTRASTS, PRIMARY_METRICS, keys)}):
        raise ValueError("contrast plotting CSV requires all registered primary outcomes and coordinates")
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
    positions = range(len(CONTRASTS))
    for metric, axis in zip(PRIMARY_METRICS, axes, strict=True):
        groups = [[replicates[contrast, metric, *key] for key in keys] for contrast in CONTRASTS]
        for position, label, values in zip(positions, CONTRASTS, groups, strict=True):
            mean = summary[label, metric]
            if not mean_matches(mean, values):
                raise ValueError(f"coordinate contrasts disagree with aggregate for {label}, {metric}")
            plotting.replicate_bar(axis, position, values, hue=220 if position < 4 else 300, depth=2, width=.68)
            if math.isnan(mean):
                axis.text(position, .02, "NA", ha="center", fontsize=7, transform=axis.get_xaxis_transform())
        axis.set_title(metric, fontweight="bold")
        axis.set_xticks(positions, tuple(CONTRASTS))
        axis.set_xlim(-.5, len(CONTRASTS) - .5)
        axis.axhline(0, color="#666666", linewidth=.8, zorder=1.5)
        plotting.focus_y(axis, groups, reference=0)
        axis.grid(axis="y", alpha=.25, zorder=1)
        axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle("Primary contrasts · Last two bars are descriptive interactions")
    figure.supylabel("Signed difference")
    plotting.save_figure(figure, output / CONTRAST_NAME, dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true", help="Render the bounded smoke coordinates")
    args = parser.parse_args()
    render_metrics(args.data, args.output, args.smoke)
    render_curves(args.data, args.output, args.smoke)
    render_contrasts(args.data, args.output, args.smoke)


if __name__ == "__main__":
    main()
