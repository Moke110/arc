"""Render the canonical Stage 2 comparison figures from plotted CSVs."""

import csv
import math
from collections import defaultdict
from pathlib import Path

from .plotting import focus_y, learning_curve, replicate_bar, save_figure

METRIC_PANELS = (
    ("Δ Spearman", "ΔSpearman"), ("wJAC-2000", "wJAC-2000"), ("JAC_top200", "JAC-top200"),
    ("Δ PCC", "ΔPCC"), ("Δ RMSE", "ΔRMSE"), ("Δ R²", "ΔR²"),
    ("Δ E-distance", "ΔE-distance"), ("Δ MMD", "ΔMMD"),
)
METRIC_NAMES = tuple(metric for metric, _ in METRIC_PANELS)
FIGURES = ("learning-curve.png", "metrics.png", "metrics-paired.png")
CURVE_FIELDS = ("arm", "fold", "seed", "epoch", "train", "validation")
METRIC_FIELDS = ("arm", "fold", "seed", "metric", "value")
PAIRED_FIELDS = (*METRIC_FIELDS[:-1], "reference", "arm_value", "reference_value", "value")


def _read(path, fields):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = csv.DictReader(stream)
        if tuple(rows.fieldnames or ()) != fields:
            raise ValueError(f"unexpected columns in {path}: {rows.fieldnames}")
        return list(rows)


def _coordinate(row):
    return row["arm"], int(row["fold"]), int(row["seed"])


def _score(value):
    if value == "":
        return None
    score = float(value)
    if not math.isfinite(score):
        raise ValueError("record undefined metrics as empty CSV fields")
    return score


def _read_metrics(data):
    scores = {}
    panels = defaultdict(set)
    for row in _read(data / "metrics.csv", METRIC_FIELDS):
        coordinate = _coordinate(row)
        metric = row["metric"]
        key = (*coordinate, metric)
        if key in scores or metric not in METRIC_NAMES:
            raise ValueError(f"duplicate or unknown metric for {key}")
        scores[key] = _score(row["value"])
        panels[coordinate].add(metric)
    if not panels or any(values != set(METRIC_NAMES) for values in panels.values()):
        raise ValueError("every fit needs all eight metrics")
    return scores, set(panels)


def _pair(scores, reference):
    if not any(arm == reference for arm, *_ in scores):
        raise ValueError(f"reference arm {reference!r} has no metrics")
    pairs = []
    for (arm, fold, seed, metric), value in scores.items():
        if arm == reference:
            continue
        baseline = scores.get((reference, fold, seed, metric))
        difference = value - baseline if value is not None and baseline is not None else None
        pairs.append((arm, fold, seed, metric, reference, value, baseline, difference))
    return pairs


def write_paired_metrics(data, reference):
    """Write the paired CSV from raw scores using exact fold and seed keys."""
    data = Path(data)
    scores, _ = _read_metrics(data)
    with (data / "metrics-paired.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(PAIRED_FIELDS)
        writer.writerows(_pair(scores, reference))


def _read_curves(data):
    curves = defaultdict(dict)
    for row in _read(data / "learning-curve.csv", CURVE_FIELDS):
        key = _coordinate(row)
        epoch = int(row["epoch"])
        losses = float(row["train"]), float(row["validation"])
        if epoch in curves[key] or not all(map(math.isfinite, losses)):
            raise ValueError(f"duplicate epoch or undefined loss for {key}, epoch {epoch}")
        curves[key][epoch] = losses
    for key, trace in curves.items():
        if sorted(trace) != list(range(1, len(trace) + 1)):
            raise ValueError(f"incomplete epoch trace for {key}")
    return curves


def _read_pairs(data, scores, reference):
    pairs = [
        (*_coordinate(row), row["metric"], row["reference"],
         _score(row["arm_value"]), _score(row["reference_value"]), _score(row["value"]))
        for row in _read(data / "metrics-paired.csv", PAIRED_FIELDS)
    ]
    if pairs != _pair(scores, reference):
        raise ValueError("paired CSV differs from the raw fold–seed scores")
    return {row[:4]: row[-1] for row in pairs}


def _format_panel(axis):
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#eeeeee")
    axis.set_axisbelow(True)


def _curve_figure(plt, curves, arms, title):
    count = len(arms)
    ncols = min(3, count)
    nrows = math.ceil(count / ncols)
    figure, axes = plt.subplots(nrows, ncols, squeeze=False,
                                figsize=(4.4 * ncols, 2.8 * nrows), layout="constrained")
    for axis, (arm, (label, _, _)) in zip(axes.flat, arms.items()):
        runs = []
        for (name, _, _), trace in sorted(curves.items()):
            if name == arm:
                epochs = sorted(trace)
                train, validation = zip(*(trace[epoch] for epoch in epochs), strict=True)
                runs.append((epochs, train, validation))
        learning_curve(axis, runs)
        axis.set_title(label)
        axis.set_xlabel("Epoch")
        axis.set_ylabel("Objective")
        axis.legend(loc="upper right", fontsize=7, frameon=False)
        _format_panel(axis)
    for axis in list(axes.flat)[count:]:
        axis.remove()
    figure.suptitle(f"{title} · learning curves")
    return figure


def _metric_figure(plt, scores, arms, title):
    figure, axes = plt.subplots(3, 3, figsize=(13.2, 10.4), layout="constrained")
    labels = [style[0] for style in arms.values()]
    for axis, (metric, label) in zip(axes.flat, METRIC_PANELS):
        groups = []
        for x, (arm, (_, hue, depth)) in enumerate(arms.items()):
            values = [value for (name, _, _, score), value in scores.items()
                      if name == arm and score == metric and value is not None]
            groups.append(values)
            if values:
                replicate_bar(axis, x, values, hue=hue, depth=depth, width=.66)
            else:
                axis.text(x, .04, "NA", ha="center", transform=axis.get_xaxis_transform())
        axis.set_xticks(range(len(arms)), labels, fontsize=8, rotation=30,
                        ha="right", rotation_mode="anchor")
        axis.set_xlim(-.5, len(arms) - .5)
        axis.axhline(0, color="#777777", lw=.6)
        focus_y(axis, groups, reference=0)
        _format_panel(axis)
        axis.set_title(label)
    axes[2, 2].remove()
    figure.suptitle(title)
    return figure


def render_comparison(data, output, *, arms, reference, title):
    """Verify the three plotted CSVs and render with the study's arm labels and colors."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data, output = Path(data), Path(output)
    scores, coordinates = _read_metrics(data)
    curves = _read_curves(data)
    if set(curves) != coordinates or set(arms) != {arm for arm, _, _ in coordinates}:
        raise ValueError("score, epoch trace and renderer arm rosters differ")
    paired = _read_pairs(data, scores, reference)
    comparison = {arm: style for arm, style in arms.items() if arm != reference}
    reference_label = arms[reference][0]
    figures = (
        _curve_figure(plt, curves, arms, title),
        _metric_figure(plt, scores, arms, f"{title} · held-out scores · mean ± SD"),
        _metric_figure(plt, paired, comparison,
                       f"{title} · paired differences vs {reference_label} · mean ± SD"),
    )
    for name, figure in zip(FIGURES, figures, strict=True):
        save_figure(figure, output / name, dpi=180)
        plt.close(figure)
