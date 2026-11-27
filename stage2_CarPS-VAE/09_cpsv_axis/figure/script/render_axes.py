"""Render the seven Step 09 figures using only their plotted CSVs."""

import argparse
from itertools import product
import math
from pathlib import Path
from statistics import stdev
import sys

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / "execution"))
from axis_protocol import AXES, DISEASE_OVERLAP, SEEDS, coordinates
from stage2_util import palette, plotting

SCORE_METRICS = ("JAC_top200", "JAC_pos100", "JAC_neg100", "Δ PCC", "Δ Spearman",
                 "Δ RMSE", "Δ E-distance", "Δ MMD", "Δ R²")
COVERAGE = "direction-correct coverage"
PLOTS = (
    ("axis-metrics", "Held-out full-axis Δ scores", False, False),
    ("axis-metrics-paired", "Paired full-axis scores versus Denoise-DEG 5k", True, False),
    ("disease-metrics", "Held-out disease-restricted Δ scores", False, True),
    ("disease-metrics-paired", "Paired disease scores versus Denoise-DEG 5k", True, True),
    ("disease-correct-coverage", "Disease direction-correct coverage", False, True),
    ("disease-correct-coverage-paired", "Paired disease coverage versus Denoise-DEG 5k", True, True),
)
FIGURES = ("learning-curves", *(spec[0] for spec in PLOTS))


def label(axis):
    family, width = axis.split("-")
    return f"{'HVG' if family == 'hvg' else 'Denoise'} {width}"


def float_or_nan(value):
    return float(value) if value else math.nan


def finite_approx(left, right):
    return (math.isnan(left) and math.isnan(right)) or math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-14)


def validate_coordinates(keys):
    if set(keys) not in (set(coordinates()), set(coordinates(True))):
        raise ValueError("plotted CSV does not cover the formal or bounded Step 09 coordinate roster")


def read_score_panels(path, metrics, paired, disease):
    groups, keys, seen = {}, set(), set()
    for row in plotting.read_rows(path):
        key = row["axis"], int(row["fold"]), int(row["seed"])
        axis = key[0]
        metric = row["metric"]
        if ((*key, metric) in seen or metric not in metrics or row["scope"] != ("disease" if disease else "axis")
                or int(row["overlap"]) != (DISEASE_OVERLAP[axis] if disease else int(axis[-2]) * 1000)):
            raise ValueError(f"invalid plotted score identity for {key}, {metric}")
        expected_zero = 0.0 if paired or metric == "Δ R²" else math.nan
        if not finite_approx(float_or_nan(row["zero_reference"]), expected_zero):
            raise ValueError(f"invalid zero reference for {key}, {metric}")
        seen.add((*key, metric))
        keys.add(key)
        groups.setdefault((axis, metric), {})[key[1:]] = tuple(
            float_or_nan(row[field]) for field in ("plotted_value", "plotted_mean", "plotted_sd"))
    validate_coordinates(keys)
    if seen != {(*key, metric) for key in keys for metric in metrics}:
        raise ValueError("plotted scores are incomplete")
    panels = {}
    for (axis, metric), group in groups.items():
        ordered = [group[pair] for pair in sorted(group)]
        values = [value for value, _, _ in ordered if math.isfinite(value)]
        mean = math.fsum(values) / len(values) if values else math.nan
        sd = stdev(values) if len(values) > 1 else math.nan
        for value, expected_mean, expected_sd in ordered:
            if math.isinf(value) or not finite_approx(mean, expected_mean) or not finite_approx(sd, expected_sd):
                raise ValueError(f"plotted values and sample SD disagree for {axis}, {metric}")
        panels[axis, metric] = [value for value, _, _ in ordered]
    return panels


def score_figure(data, output, name, title, paired, disease):
    coverage = "coverage" in name
    metrics = (COVERAGE,) if coverage else SCORE_METRICS
    panels = read_score_panels(data / f"{name}.csv", metrics, paired, disease)
    dimensions = ((1, 1), (11, 5)) if coverage else ((3, 3), (19, 12))
    figure, axes = plt.subplots(*dimensions[0], figsize=dimensions[1], squeeze=False, layout="constrained")
    for axis_plot, metric in zip(axes.flat, metrics, strict=True):
        groups = []
        for i, axis in enumerate(AXES):
            values = panels[axis, metric]
            groups.append(values)
            hue = palette.HUES[0 if axis.startswith("hvg") else 1]
            depth = int(axis[-2]) - 2
            plotting.replicate_bar(axis_plot, i, values, hue=hue, depth=depth, width=.75)
        zero = 0 if paired or metric == "Δ R²" else None
        if zero is not None:
            axis_plot.axhline(zero, color=".32", lw=.8, zorder=1)
        axis_plot.set_xticks(range(len(AXES)), [
            label(axis) + (f"\n{DISEASE_OVERLAP[axis]:,} genes" if disease else "") for axis in AXES
        ], fontsize=8, rotation=24, ha="right")
        axis_plot.set_title(metric.replace("JAC_", "JAC ") if not coverage else (
            "Difference (percentage points)" if paired else "Coverage of 3,481 genes (%)"))
        plotting.focus_y(axis_plot, groups, reference=zero)
        axis_plot.grid(axis="y", alpha=.2)
    figure.suptitle(title + ("\nY-axis limits may omit zero." if not paired else ""), fontsize=15)
    figure.savefig(output / f"{name}.png", dpi=140)
    plt.close(figure)


def learning_figure(data, output):
    curves = {}
    for row in plotting.read_rows(data / "learning-curves.csv"):
        key = row["axis"], int(row["fold"]), int(row["seed"])
        epoch = int(row["epoch"])
        curve = curves.setdefault(key, {})
        losses = float(row["train_loss"]), float(row["validation_loss"])
        if epoch in curve or not all(map(math.isfinite, losses)):
            raise ValueError(f"invalid objective trace for {key}")
        curve[epoch] = losses
    validate_coordinates(curves)
    figure, axes = plt.subplots(2, 4, figsize=(16, 8), squeeze=False, layout="constrained")
    for panel, axis in zip(axes.flat, AXES, strict=True):
        for fold, seed in product(range(3), SEEDS):
            curve = curves.get((axis, fold, seed))
            if curve is None:
                continue
            epochs = sorted(curve)
            if epochs != list(range(1, 1 + len(epochs))):
                raise ValueError("objective trace epochs are not contiguous from one")
            train, valid = zip(*(curve[epoch] for epoch in epochs), strict=True)
            style = {"hue": palette.HUES[fold], "depth": SEEDS.index(seed) + 1, "linewidth": 1.15}
            plotting.line(panel, epochs, train, linestyle="-", **style)
            plotting.line(panel, epochs, valid, linestyle="--", **style)
        panel.set_title(label(axis))
        panel.set_xlabel("Epoch")
        panel.set_ylabel("Objective per condition")
        panel.grid(alpha=.2)
    legend = [Line2D([0], [0], color="black", linestyle=style, label=label_text)
              for label_text, style in (("Training", "-"), ("Validation", "--"))]
    for fold in range(3):
        color = palette.color_pair(palette.HUES[fold], 2).main
        legend.append(Line2D([0], [0], color=color, label=f"Fold {fold}"))
    for depth, seed in enumerate(SEEDS, start=1):
        color = palette.color_pair(palette.HUES[0], depth).main
        legend.append(Line2D([0], [0], color=color, label=f"Seed {seed} shade"))
    figure.legend(handles=legend, loc="outside lower center", ncols=len(legend), fontsize=9)
    figure.suptitle("CarPS-VAE UniMol objective by prediction axis", fontsize=15)
    figure.savefig(output / "learning-curves.png", dpi=140)
    plt.close(figure)


def render(data, output):
    output.mkdir(parents=True, exist_ok=True)
    learning_figure(data, output)
    for spec in PLOTS:
        score_figure(data, output, *spec)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == "__main__":
    main()
