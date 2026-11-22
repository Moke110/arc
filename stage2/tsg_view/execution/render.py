#!/usr/bin/env python3
"""Render the `tsg_view` figures from ``figure/data/*.csv`` alone.

The renderer is deliberately self-contained: it imports only the standard library and
matplotlib, reads only the CSVs named below from ``--data``, validates their complete
grid, and draws the values it is given.  It never opens a run directory, a matrix, a
checkpoint or the truth product, and it recomputes no scientific metric (means, spreads
and interval bounds are stored in the CSVs by :mod:`report`).

Byte-identical re-render is enforced by construction: Agg backend fixed before pyplot,
explicit figure size and subplot spacing, fixed dpi, ``facecolor='white'`` and
``metadata={'Software': None}`` (no matplotlib version chunk, no timestamp).  The
companion test renders twice in fresh ``python -I`` processes with isolated
``MPLCONFIGDIR`` into different directories and compares PNG sha256.

Suites keep the two-pass publication honest: ``--suite tenv`` writes only the tenv
panels and records Tahoe as explicitly pending; ``--suite tahoe`` writes only the
Tahoe panels; ``--suite cross`` writes only ``views_effect.png``; ``--suite all``
requires both sources.  Nothing in one suite deletes or rewrites another suite's PNGs.

CLI::

    render.py --data <figure/data> --output <directory> --suite tenv|tahoe|all|cross
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

SCHEMA = "tsg-view-render-v1"
RENDERER_VERSION = "tsg-view-render-v1"
DPI = 150
SOURCES = ("tenv", "tahoe_bulk")
VIEW_ORDER = {
    "tenv": ("full", "merge", "solo", "full-pooled-label", "solo-deepest"),
    "tahoe_bulk": ("full", "merge", "solo"),
}
SEEDS = {"tenv": {view: (42,) if view == "solo-deepest" else (42, 43)
                  for view in VIEW_ORDER["tenv"]},
         "tahoe_bulk": {view: (42, 43) for view in VIEW_ORDER["tahoe_bulk"]}}
FOLDS = (0, 1, 2)
METRICS = (("delta_pcc", "\u0394 Pearson r (PCC)"),
           ("delta_spearman", "\u0394 Spearman r"),
           ("delta_rmse", "\u0394 RMSE"),
           ("delta_r2", "\u0394 R\u00b2"))
METRIC_LABELS = dict(METRICS)
METRIC_ORDER = tuple(metric for metric, _ in METRICS)
COLORS = ("#1f4e79", "#c55a11", "#2f7d32", "#7b3f9d", "#b8860b", "#3a7ca5")
SOURCE_NAMES = {"tenv": "tenv", "tahoe_bulk": "Tahoe"}

FILES = {
    "learning_curves": ("source", "view", "arm_type", "fold", "seed", "epoch", "train_loss",
                        "validation_mse", "best_validation"),
    "learning_curves_mean": ("source", "view", "arm_type", "epoch", "mean_train_loss",
                             "mean_validation_mse", "n_runs"),
    "metrics_bars": ("source", "view", "arm_type", "metric", "headline", "fold_sd", "seed_sd",
                     "n_folds", "n_seeds", "sd_estimand", "mean_cell", "mean_cell_fold_sd"),
    "metrics_points": ("source", "view", "arm_type", "fold", "seed", "metric", "value",
                       "x_offset"),
    "metrics_baseline": ("source", "metric", "mean_cell", "mean_cell_fold_sd", "n_folds",
                         "provenance_sha256"),
    "metrics_per_dose_bars": ("source", "fold", "held_dose_um", "held_dose_label", "view",
                              "arm_type", "metric", "value", "seed_sd", "n_seeds", "n_samples",
                              "dose_fold_confounded", "design_note"),
    "metrics_per_dose_points": ("source", "fold", "held_dose_um", "view", "arm_type", "metric",
                                "seed", "value", "x_offset"),
    "views_effect": ("source", "view", "arm_type", "metric", "delta_mean", "low", "high",
                     "fold_sd", "seed_sd", "n_pairs", "status", "interval_estimator"),
}


class RenderError(ValueError):
    """The figure data are not the complete grid the figure contract promises."""


# --------------------------------------------------------------------------- #
# reading and validation
# --------------------------------------------------------------------------- #
def read_rows(data, name):
    path = Path(data) / f"{name}.csv"
    if not path.is_file():
        raise RenderError(f"missing figure data: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != FILES[name]:
            raise RenderError(f"{name}.csv header does not match the renderer contract")
        rows = list(reader)
    if not rows:
        raise RenderError(f"{name}.csv carries no rows")
    return rows


def number(row, name):
    try:
        return float(row[name])
    except (KeyError, TypeError, ValueError) as error:
        raise RenderError(f"non-numeric {name!r} in figure data") from error


def optional_number(row, name):
    value = row.get(name, "")
    if value is None or str(value).strip() == "":
        return None
    return number(row, name)


def integer(row, name):
    value = number(row, name)
    if value != int(value):
        raise RenderError(f"non-integer {name!r} in figure data")
    return int(value)


def checked_index(rows, keys, expected, name):
    index = {}
    for row in rows:
        key = tuple(row[column] for column in keys)
        if key in index:
            raise RenderError(f"duplicate {keys} group in {name}: {key}")
        index[key] = row
    missing = sorted(set(expected) - set(index))
    extra = sorted(set(index) - set(expected))
    if missing or extra:
        raise RenderError(f"{name} grid mismatch: missing={missing[:4]} extra={extra[:4]}")
    return index


def source_rows(rows, source):
    return [row for row in rows if row["source"] == source]


def available_sources(data, name="metrics_bars"):
    rows = read_rows(data, name)
    return [source for source in SOURCES if any(row["source"] == source for row in rows)]


# --------------------------------------------------------------------------- #
# drawing primitives
# --------------------------------------------------------------------------- #
def save(figure, output, name):
    path = Path(output) / f"{name}.png"
    figure.savefig(path, dpi=DPI, facecolor="white", metadata={"Software": None})
    plt.close(figure)
    return path


def style(axis):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.grid(axis="y", color="0.9", linewidth=0.6)
    axis.set_axisbelow(True)


def _subplots(rows, columns, width, height, **adjust):
    figure, axes = plt.subplots(rows, columns, figsize=(width, height), squeeze=False)
    figure.subplots_adjust(**adjust)
    return figure, axes


# --------------------------------------------------------------------------- #
# learning curves
# --------------------------------------------------------------------------- #
def render_learning_curves(data, source, output):
    rows = source_rows(read_rows(data, "learning_curves"), source)
    means = source_rows(read_rows(data, "learning_curves_mean"), source)
    views = VIEW_ORDER[source]
    expected = {(view, str(fold), str(seed))
                for view in views for fold in FOLDS for seed in SEEDS[source][view]}
    actual = {(row["view"], row["fold"], row["seed"]) for row in rows}
    missing = sorted(expected - actual)
    extra = sorted(actual - expected)
    if missing or extra:
        raise RenderError(f"learning_curves grid mismatch: missing={missing[:4]} "
                          f"extra={extra[:4]}")
    curves = {}
    for key in sorted(actual):
        view, fold, seed = key
        curve = sorted((row for row in rows
                        if row["view"] == view and row["fold"] == fold and row["seed"] == seed),
                       key=lambda row: integer(row, "epoch"))
        epochs = [integer(row, "epoch") for row in curve]
        if epochs != list(range(1, len(epochs) + 1)):
            raise RenderError(f"learning curve is not contiguous: {view}/{fold}/{seed}")
        best = [row for row in curve if integer(row, "best_validation") == 1]
        if len(best) != 1:
            raise RenderError(f"learning curve has {len(best)} best epochs: {view}/{fold}/{seed}")
        curves[(view, int(fold), int(seed))] = curve
    mean_curves = {}
    for view in views:
        block = sorted((row for row in means if row["view"] == view),
                       key=lambda row: integer(row, "epoch"))
        epochs = [integer(row, "epoch") for row in block]
        if epochs != list(range(1, len(epochs) + 1)):
            raise RenderError(f"mean learning curve is not contiguous: {view}")
        mean_curves[view] = block
    maximum = max(integer(row, "epoch") for row in rows)
    figure, axes = _subplots(1, len(views), 3.4 * len(views), 3.2,
                             left=0.05, right=0.99, top=0.86, bottom=0.13, wspace=0.24)
    for column, view in enumerate(views):
        axis = axes[0][column]
        style(axis)
        for fold in FOLDS:
            for seed in SEEDS[source][view]:
                curve = curves[(view, fold, int(seed))]
                axis.plot([integer(row, "epoch") for row in curve],
                          [number(row, "validation_mse") for row in curve],
                          color=COLORS[fold % len(COLORS)], linewidth=0.8, alpha=0.45)
                best = [row for row in curve if integer(row, "best_validation") == 1][0]
                axis.plot([integer(best, "epoch")], [number(best, "validation_mse")],
                          marker="o", markersize=3.2,
                          color=COLORS[fold % len(COLORS)], linestyle="none")
        block = mean_curves[view]
        if block:
            axis.plot([integer(row, "epoch") for row in block],
                      [number(row, "mean_validation_mse") for row in block],
                      color="black", linewidth=1.6)
        axis.set_title(view, fontsize=10)
        axis.set_xlabel("epoch", fontsize=9)
        axis.set_xticks([1, maximum])
        if column == 0:
            axis.set_ylabel("validation MSE", fontsize=9)
        if column == 0:
            handles = [Line2D([0], [0], color=COLORS[fold % len(COLORS)], linewidth=1.0,
                              label=f"fold {fold}") for fold in FOLDS]
            handles.append(Line2D([0], [0], color="black", linewidth=1.6, label="run mean"))
            axis.legend(handles=handles, loc="upper right", fontsize=7, frameon=False)
    figure.suptitle(f"{SOURCE_NAMES[source]} learning curves (validation MSE; "
                    "dots mark each run's best epoch)", fontsize=11)
    return save(figure, output, f"learning_curves_{source}")


# --------------------------------------------------------------------------- #
# metrics bars
# --------------------------------------------------------------------------- #
def _validate_metric_grid(rows, source, name):
    views = VIEW_ORDER[source]
    expected = {(view, metric) for view in views for metric in METRIC_ORDER}
    return checked_index(rows, ("view", "metric"), expected, name)


def render_metrics(data, source, output):
    bars = source_rows(read_rows(data, "metrics_bars"), source)
    points = source_rows(read_rows(data, "metrics_points"), source)
    baseline = source_rows(read_rows(data, "metrics_baseline"), source)
    views = VIEW_ORDER[source]
    bar_index = _validate_metric_grid(bars, source, "metrics_bars")
    expected_points = {(view, str(fold), str(seed), metric)
                       for view in views for fold in FOLDS for seed in SEEDS[source][view]
                       for metric in METRIC_ORDER}
    checked_index(points, ("view", "fold", "seed", "metric"), expected_points, "metrics_points")
    expected_baseline = {(metric,) for metric in METRIC_ORDER}
    checked_index(baseline, ("metric",), expected_baseline, "metrics_baseline")
    baseline_value = {row["metric"]: number(row, "mean_cell") for row in baseline}
    figure, axes = _subplots(1, len(METRICS), 4.0 * len(METRICS), 3.9,
                             left=0.05, right=0.99, top=0.82, bottom=0.22, wspace=0.20)
    positions = list(range(len(views)))
    for column, (metric, label) in enumerate(METRICS):
        axis = axes[0][column]
        style(axis)
        for position, view in enumerate(views):
            row = bar_index[(view, metric)]
            headline = number(row, "headline")
            fold_sd = optional_number(row, "fold_sd")
            seed_sd = optional_number(row, "seed_sd")
            axis.bar([position], [headline], width=0.46, color=COLORS[position % len(COLORS)],
                     alpha=0.75,
                     hatch="///" if row["arm_type"].startswith("solo-deepest") else None,
                     edgecolor="black", linewidth=0.5)
            if fold_sd is not None:
                axis.errorbar([position - 0.11], [headline], yerr=[fold_sd], fmt="none",
                              ecolor="black", elinewidth=1.2, capsize=2.5, capthick=1.2,
                              label="fold SD" if position == 0 and column == 0 else None)
            if seed_sd is not None:
                axis.errorbar([position + 0.11], [headline], yerr=[seed_sd], fmt="none",
                              ecolor="0.45", elinewidth=1.2, capsize=2.5, capthick=1.2,
                              linestyle="--",
                              label="seed SD" if position == 0 and column == 0 else None)
            run_points = [row for row in points
                          if row["view"] == view and row["metric"] == metric]
            for point in run_points:
                axis.scatter(number(point, "x_offset") + position, number(point, "value"),
                             s=9, color=COLORS[integer(point, "fold") % len(COLORS)],
                             alpha=0.75, zorder=3)
            spread = max([value for value in (fold_sd, seed_sd) if value is not None] or [0.0])
            suffix = "" if int(row["n_seeds"]) > 1 else " (1 seed)"
            axis.text(position, headline + spread * 0.25, f"{headline:.3f}{suffix}",
                      ha="center", va="bottom", fontsize=5.5, rotation=90)
        axis.axhline(baseline_value[metric], color="0.35", linestyle="--", linewidth=1.0,
                     label="Mean_cell" if column == 0 else None)
        axis.axhline(0, color="0.8", linewidth=0.7)
        axis.set_xticks(positions)
        axis.set_xticklabels(views, rotation=35, ha="right", fontsize=8)
        axis.set_title(label, fontsize=10)
        if column == 0:
            axis.set_ylabel("equal-fold mean", fontsize=9)
            axis.legend(loc="best", fontsize=7, frameon=False)
    figure.suptitle(f"{SOURCE_NAMES[source]} view scores: equal-fold mean with fold SD "
                    "(black) and seed SD (grey) drawn separately; dots are runs; "
                    "dashed line is Mean_cell", fontsize=10)
    return save(figure, output, f"metrics_{source}")


# --------------------------------------------------------------------------- #
# per-dose panels
# --------------------------------------------------------------------------- #
def _validate_doses(rows, source):
    views = VIEW_ORDER[source]
    folds = {}
    for fold in FOLDS:
        block = [row for row in rows if integer(row, "fold") == fold]
        doses = {number(row, "held_dose_um") for row in block}
        if len(doses) != 1:
            raise RenderError(f"fold {fold} does not hold exactly one dose: {sorted(doses)}")
        for row in block:
            if row["dose_fold_confounded"].lower() != "true":
                raise RenderError("per-dose rows must state that dose and fold are confounded")
        folds[fold] = doses.pop()
    if len(set(folds.values())) != len(FOLDS):
        raise RenderError("per-dose panels do not cover three distinct held doses")
    expected = {(str(fold), view, metric) for fold in FOLDS for view in views
                for metric in METRIC_ORDER}
    checked_index(rows, ("fold", "view", "metric"), expected, "metrics_per_dose_bars")
    return folds


def render_metrics_per_dose(data, source, output):
    bars = source_rows(read_rows(data, "metrics_per_dose_bars"), source)
    points = source_rows(read_rows(data, "metrics_per_dose_points"), source)
    views = VIEW_ORDER[source]
    doses = _validate_doses(bars, source)
    expected_points = {(str(fold), view, metric, str(seed))
                       for fold in FOLDS for view in views for metric in METRIC_ORDER
                       for seed in SEEDS[source][view]}
    checked_index(points, ("fold", "view", "metric", "seed"), expected_points,
                  "metrics_per_dose_points")
    sample_counts = {}
    for fold in FOLDS:
        values = {integer(row, "n_samples") for row in bars if integer(row, "fold") == fold}
        if len(values) != 1:
            raise RenderError(f"fold {fold} carries inconsistent per-dose sample counts")
        sample_counts[fold] = values.pop()
    figure, axes = _subplots(len(METRICS), len(FOLDS), 4.3 * len(FOLDS), 2.9 * len(METRICS),
                             left=0.06, right=0.985, top=0.94, bottom=0.10, wspace=0.28,
                             hspace=0.48)
    positions = list(range(len(views)))
    for row_index, (metric, label) in enumerate(METRICS):
        for column, fold in enumerate(FOLDS):
            axis = axes[row_index][column]
            style(axis)
            for position, view in enumerate(views):
                row = [candidate for candidate in bars
                       if integer(candidate, "fold") == fold
                       and candidate["view"] == view and candidate["metric"] == metric][0]
                value = number(row, "value")
                axis.bar([position], [value], width=0.52,
                         color=COLORS[position % len(COLORS)], alpha=0.75,
                         hatch="///" if row["arm_type"].startswith("solo-deepest") else None,
                         edgecolor="black", linewidth=0.5)
                seed_sd = optional_number(row, "seed_sd")
                if seed_sd is not None:
                    axis.errorbar([position], [value], yerr=[seed_sd], fmt="none",
                                  ecolor="0.35", elinewidth=1.0, capsize=2.0, capthick=1.0)
                fold_points = [point for point in points
                               if integer(point, "fold") == fold and point["view"] == view
                               and point["metric"] == metric]
                for point in fold_points:
                    axis.scatter(number(point, "x_offset") + position, number(point, "value"),
                                 s=9, color="black", alpha=0.7, zorder=3)
            axis.axhline(0, color="0.8", linewidth=0.7)
            axis.set_xticks(positions)
            axis.set_xticklabels(views, rotation=35, ha="right", fontsize=7)
            axis.set_title(f"fold {fold}: {doses[fold]:g} \u00b5M \u00b7 "
                           f"{sample_counts[fold]} test profiles", fontsize=8.5)
            if column == 0:
                axis.set_ylabel(label, fontsize=8.5)
    figure.legend(handles=[Line2D([0], [0], marker="o", color="black", linestyle="none",
                                  markersize=3.5, label="run (dot)"),
                           Line2D([0], [0], color="0.35", linewidth=1.0,
                                  label="seed SD (whisker)")],
                  loc="lower center", ncol=2, fontsize=7, frameon=False)
    figure.suptitle(f"{SOURCE_NAMES[source]} per-dose panels: dose and fold are confounded by "
                    "design (each fold holds exactly one dose); bars are seed means, dots are "
                    "runs", fontsize=10)
    return save(figure, output, f"metrics_per_dose_{source}")


# --------------------------------------------------------------------------- #
# cross-source effect
# --------------------------------------------------------------------------- #
def render_views_effect(data, output):
    rows = read_rows(data, "views_effect")
    sources = [source for source in SOURCES if any(row["source"] == source for row in rows)]
    if set(sources) != set(SOURCES):
        raise RenderError(
            "cross-source views_effect requires both sources; pending: "
            f"{sorted(set(SOURCES) - set(sources))}")
    views = [view for view in VIEW_ORDER["tenv"] if view != "full"]
    index = {(row["source"], row["view"], row["metric"]): row for row in rows}
    figure, axes = _subplots(1, len(METRICS), 4.0 * len(METRICS), 4.0,
                             left=0.05, right=0.99, top=0.80, bottom=0.24, wspace=0.20)
    positions = list(range(len(views)))
    for column, (metric, label) in enumerate(METRICS):
        axis = axes[0][column]
        style(axis)
        for position, view in enumerate(views):
            for offset, source in ((-0.18, "tenv"), (0.18, "tahoe_bulk")):
                row = index.get((source, view, metric))
                if row is None:
                    axis.text(position + offset, 0, "not fitted", rotation=90, fontsize=6,
                              ha="center", va="center", color="0.4")
                    continue
                if row["status"] != "ok":
                    axis.text(position + offset, 0, row["status"], rotation=90, fontsize=6,
                              ha="center", va="center", color="0.4")
                    continue
                mean = number(row, "delta_mean")
                low, high = number(row, "low"), number(row, "high")
                axis.bar([position + offset], [mean], width=0.30,
                         color=COLORS[0] if source == "tenv" else COLORS[1], alpha=0.8,
                         edgecolor="black", linewidth=0.4,
                         hatch="///" if row["arm_type"].startswith("solo-deepest") else None)
                axis.errorbar([position + offset], [mean],
                              yerr=[[mean - low], [high - mean]], fmt="none",
                              ecolor="black", elinewidth=1.0, capsize=2.0, capthick=1.0)
        axis.axhline(0, color="0.8", linewidth=0.7)
        axis.set_xticks(positions)
        axis.set_xticklabels([view for view in views], rotation=35, ha="right", fontsize=8)
        axis.set_title(f"{label}: paired difference vs full", fontsize=9)
        if column == 0:
            axis.set_ylabel("\u0394(view) - \u0394(full), equal-fold mean", fontsize=8.5)
            axis.legend(handles=[Line2D([0], [0], color=COLORS[0], linewidth=6, label="tenv"),
                                 Line2D([0], [0], color=COLORS[1], linewidth=6, label="Tahoe")],
                        loc="best", fontsize=8, frameon=False)
    figure.suptitle("Cross-source view effect: paired view - full differences with a "
                    "fold+seed interval; solo-deepest is hatched and reported separately",
                    fontsize=10)
    return save(figure, output, "views_effect")


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #
def render(data, output, suite="all"):
    data = Path(data)
    output = Path(output)
    source_list = available_sources(data)
    if suite in ("tenv", "tahoe_bulk"):
        wanted = (suite,)
        if suite not in source_list:
            raise RenderError(f"--suite {suite} requested but the figure data carry no "
                              f"{suite} rows (available: {source_list})")
    elif suite == "all":
        wanted = SOURCES
        missing = [source for source in SOURCES if source not in source_list]
        if missing:
            raise RenderError(f"--suite all requires both sources; pending: {missing}")
    elif suite == "cross":
        wanted = ()
    else:
        raise RenderError(f"unknown suite: {suite}")
    output.mkdir(parents=True, exist_ok=True)
    figures = {}
    if suite == "cross":
        for name in ("views_effect.png",):
            figures[f"views_effect"] = render_views_effect(data, output)
    else:
        for source in wanted:
            figures[f"learning_curves_{source}"] = render_learning_curves(data, source, output)
            figures[f"metrics_{source}"] = render_metrics(data, source, output)
            figures[f"metrics_per_dose_{source}"] = render_metrics_per_dose(data, source, output)
        if suite == "all":
            figures["views_effect"] = render_views_effect(data, output)
    pending = {}
    if suite == "tenv":
        pending["tahoe_bulk"] = "pending: Tahoe arms are absent from the report data (tenv suite first)"
    elif suite == "tahoe_bulk":
        pending["tenv"] = "not part of this suite; accepted tenv panels are left untouched"
    return figures, pending


def status_value(data, suite, figures, pending):
    csv_inputs = {}
    for name in FILES:
        path = Path(data) / f"{name}.csv"
        if path.is_file():
            csv_inputs[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "schema": SCHEMA,
        "renderer_version": RENDERER_VERSION,
        "suite": suite,
        "figures": {f"{name}.png": hashlib.sha256(path.read_bytes()).hexdigest()
                    for name, path in sorted(figures.items())},
        "csv_inputs": dict(sorted(csv_inputs.items())),
        "pending": dict(sorted(pending.items())),
        "matplotlib": matplotlib.__version__,
        "renderer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "dpi": DPI,
        "metadata_software": None,
        "estimands": {
            "metrics": "bar = equal-fold mean; fold SD and seed SD are two separate error bars",
            "per_dose": "each fold holds exactly one dose; dose and fold are confounded by design",
            "views_effect": "paired view - full difference with a fold+seed interval; "
                            "solo-deepest is reported separately",
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).resolve().parents[1] / "figure" / "data")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--suite", choices=("tenv", "tahoe", "tahoe_bulk", "all", "cross"),
                        default="all")
    args = parser.parse_args(argv)
    suite = "tahoe_bulk" if args.suite == "tahoe" else args.suite
    figures, pending = render(args.data, args.output, suite)
    value = status_value(args.data, suite, figures, pending)
    (Path(args.output) / "render_status.json").write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n")
    sys.stdout.write(json.dumps({"figures": sorted(value["figures"]),
                                 "pending": value["pending"]}, indent=2,
                                sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
