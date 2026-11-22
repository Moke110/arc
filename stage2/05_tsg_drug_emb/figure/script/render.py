"""Render the source-local cold-drug figures using only plotted CSV data."""
import argparse
import csv
import json
import math
from pathlib import Path
from runpy import run_path

FIGURE = Path(__file__).resolve().parents[1]
LEGACY = (FIGURE.parents[1] / "06_tsg_dose" /
          "execution/historical/report_axis_local/script/render_axis_local.py")
SCENARIO = "cold-drug"
LABELS = {"tenv": "tenv", "tahoe_bulk": "Tahoe"}
ARMS = ("ECFP4C", "ECFP6C", "FCFP4C", "FCFP6C", "E3FP", "KPGT", "UniMol")
METHODS = ("Mean_cell", *ARMS)
COLORS = dict(zip(METHODS, ("#474747", "#465970", "#9c6032", "#167c94", "#ad4a5c",
                           "#387048", "#705ba7", "#a45f00"), strict=True))
TABLES = ("runs", "summary", "curves", "dose")
DOSES = {"tenv": (.2, 1., 10.), "tahoe_bulk": (.05, .5, 5.)}


def legacy():
    return run_path(str(LEGACY))


def load(base, data, source):
    result = {name: base["read_rows"](data / f"{name}.csv") for name in TABLES}
    for rows in result.values():
        if any((row["source"], row["scenario"], row["seed"]) != (source, SCENARIO, "42")
               for row in rows):
            raise ValueError("unexpected plot source, scenario or seed")
    metrics = base["METRICS"]
    for table, fields, coordinates in (
            ("runs", ("method", "fold"), {str(fold) for fold in range(3)}),
            ("summary", ("method", "metric"), set(metrics)),
            ("dose", ("method", "dose_um"), {str(dose) for dose in DOSES[source]})):
        index = base["by_key"](result[table], *fields)
        if (set(index) != {(arm, coordinate) for arm in METHODS for coordinate in coordinates}
                or any(len(rows) != 1 for rows in index.values())):
            raise ValueError(f"duplicate or incomplete {table} coordinates")
    curves = base["by_key"](result["curves"], "method", "fold")
    if (set(curves) != {(arm, str(fold)) for arm in ARMS for fold in range(3)}
            or any(sorted(int(row["epoch"]) for row in rows) != list(range(1, len(rows) + 1))
                   for rows in curves.values())):
        raise ValueError("duplicate or incomplete curve coordinates")
    return result


def prepare_axes(data, source="tenv"):
    base = legacy()
    tables = load(base, data, source)
    export = base["export_axes"]
    export.__globals__.update(SOURCES=(source,), SCENARIOS=(SCENARIO,), METHODS=ARMS)
    axes = export({f"round1_{name}.csv": rows for name, rows in tables.items()})
    # Preserve the published dose report's zero-based aggregate bar axes.
    for row in axes:
        if row["figure"] != "overall" or row["panel"] == "r2":
            continue
        metric = row["panel"]
        spreads = [value for value in tables["summary"] if value["metric"] == metric]
        plotted = base["overview_range"](tables["runs"], spreads, metric)
        segments = json.loads(row["segments"])
        if segments[-1]["upper"] > 0 and all(value >= 0 for value in plotted if math.isfinite(value)):
            row["segments"] = json.dumps([base["linear_segment"](
                0., segments[-1]["upper"], show_zero=True)])
    with (data / "axis_segments.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(axes[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(axes)


def overall(base, tables, axes, output, source, prefix):
    dots = base["by_key"](tables["runs"], "method")
    spreads = base["by_key"](tables["summary"], "method", "metric")
    number = base["number"]
    figure, panels = base["plt"].subplots(1, 4, figsize=(19, 4.8))
    for ax, metric, title in zip(panels, base["METRICS"], base["TITLES"], strict=True):
        scale = axes["overall", metric, "metric"]
        scale.apply(ax)
        for pos, method in enumerate(METHODS):
            (spread,) = spreads[method, metric]
            mean, sd = number(spread["mean"]), number(spread["sd"])
            color = COLORS[method]
            if math.isfinite(mean):
                scale.bar(ax, pos, mean, color)
            if math.isfinite(sd):
                lo, mid, hi = scale.map([mean - sd, mean, mean + sd])
                ax.errorbar(pos, mid, yerr=[[mid - lo], [hi - mid]],
                            color=color, capsize=3, lw=1.5)
            rows = sorted(dots[method,], key=lambda row: int(row["fold"]))
            ax.scatter([pos - .13, pos, pos + .13], scale.map([number(row[metric]) for row in rows]),
                       color=color, s=24, edgecolors="white", linewidths=.3, zorder=4)
        ax.set_title(title.rsplit(" ", 1)[0])
        ax.set_xticks(range(len(METHODS)), METHODS, rotation=32, ha="right")
        ax.set_xlim(-.55, len(METHODS) - .45)
        ax.grid(axis="y", alpha=.15)
    figure.suptitle(f"{LABELS[source]} cold-drug · equal-fold test means ± sample SD", fontsize=17)
    figure.text(.5, .035, "Each dot is one outer fold. SD shows fold dispersion for seed 42.",
                ha="center")
    figure.subplots_adjust(left=.04, right=.988, top=.82, bottom=.27, wspace=.28)
    base["save_figure"](figure, output, f"{prefix}metrics.png")


def learning(base, tables, axes, output, source, prefix):
    curves = base["by_key"](tables["curves"], "method", "fold")
    train_colors, valid_colors = base["TRAIN_COLORS"], base["VALID_COLORS"]
    number = base["number"]
    figure, panels = base["plt"].subplots(1, 8, figsize=(35, 5.1))
    overview = panels[-1]
    common = axes["learning", "all-methods", "mse"]
    common.apply(overview, color="#2e4150")
    for ax, method in zip(panels[:-1], ARMS, strict=True):
        right = ax.twinx()
        loss = axes["learning", method, "loss"]
        mse = axes["learning", method, "mse"]
        loss.apply(ax, color=train_colors[0])
        right.spines["right"].set_visible(True)
        mse.apply(right, side="right", zero_line=False, color=valid_colors[0])
        right.axhline(mse.map([0])[0], color=valid_colors[0], lw=.55, alpha=.55)
        for fold in range(3):
            rows = sorted(curves[method, str(fold)], key=lambda row: int(row["epoch"]))
            epochs = [int(row["epoch"]) for row in rows]
            loss.plot(ax, epochs, [number(row["train_loss"]) for row in rows], "--",
                      color=train_colors[fold], lw=1.)
            values = [number(row["validation_mse"]) for row in rows]
            mse.plot(right, epochs, values, color=valid_colors[fold], lw=1.)
            common.plot(overview, epochs, values, color=COLORS[method], alpha=.55, lw=.95,
                        label=method if fold == 0 else None)
        ax.set_title(method)
        ax.set_xlabel("Epoch")
        ax.set_ylabel("Native train loss", color=train_colors[0], fontsize=9)
        right.set_ylabel("Validation MSE", color=valid_colors[0], fontsize=9)
    overview.set_title("All representations")
    overview.set_xlabel("Epoch")
    overview.set_ylabel("Validation MSE · log data segment", fontsize=9)
    overview.legend(loc="upper right", fontsize=7)
    handles = [
        base["Line2D"]([0], [0], color=color, ls=style, label=f"Fold {fold} · {name}")
        for name, colors, style in (("train", train_colors, "--"), ("validation", valid_colors, "-"))
        for fold, color in enumerate(colors)
    ]
    figure.legend(handles=handles, loc="lower center", ncol=6, bbox_to_anchor=(.5, .018))
    figure.suptitle(f"{LABELS[source]} cold-drug · native training loss and validation MSE by epoch · seed 42",
                    fontsize=18)
    figure.subplots_adjust(left=.033, right=.989, top=.81, bottom=.21, wspace=.78)
    base["save_figure"](figure, output, f"{prefix}learning_curves.png")


def by_dose(base, tables, axes, output, source, prefix):
    doses = base["by_key"](tables["dose"], "method")
    grid = DOSES[source]
    number = base["number"]
    figure, panels = base["plt"].subplots(1, 4, figsize=(19, 4.8))
    for ax, metric, title in zip(panels, base["METRICS"], base["TITLES"], strict=True):
        scale = axes["dose", metric, "metric"]
        scale.apply(ax)
        for method in METHODS:
            rows = sorted(doses[method,], key=lambda row: float(row["dose_um"]))
            ax.plot(grid, scale.map([number(row[metric]) for row in rows]), "-o",
                    color=COLORS[method], label=method, lw=1.5, ms=4)
        ax.set_xscale("log")
        ax.set_xticks(grid, [f"{value:g}" for value in grid])
        ax.set_xlabel("Measured dose (µM)")
        ax.set_title(title.rsplit(" ", 1)[0])
        ax.grid(axis="y", alpha=.15)
    handles, labels = panels[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=len(labels), bbox_to_anchor=(.5, .014))
    figure.suptitle(f"{LABELS[source]} cold-drug · test performance by measured dose", fontsize=17)
    figure.text(.5, .095, "Each test sample has equal weight within its measured dose.", ha="center")
    figure.subplots_adjust(left=.04, right=.988, top=.82, bottom=.235, wspace=.28)
    base["save_figure"](figure, output, f"{prefix}metrics_cold_drug_dose.png")


def render(data, output, source="tenv", prefixed=False):
    prefix = f"{source.removesuffix('_bulk')}_" if prefixed else ""
    base = legacy()
    tables = load(base, data, source)
    axes = {}
    for row in base["read_rows"](data / "axis_segments.csv"):
        key = (row["figure"], row["panel"], row["ordinate"])
        if (row["source"], row["scenario"]) != (source, SCENARIO) or key in axes:
            raise ValueError("unexpected or duplicate axis coordinates")
        axes[key] = base["AxisScale"](json.loads(row["segments"]))
    output.mkdir(parents=True, exist_ok=True)
    with base["plt"].rc_context({"font.size": 10, "axes.spines.top": False,
                                  "axes.spines.right": False}):
        overall(base, tables, axes, output, source, prefix)
        learning(base, tables, axes, output, source, prefix)
        by_dose(base, tables, axes, output, source, prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=tuple(LABELS), default="tenv")
    parser.add_argument("--prefixed", action="store_true")
    parser.add_argument("--data", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    default_data = FIGURE / "data"
    if args.source == "tahoe_bulk":
        default_data /= "tahoe"
    render(args.data or default_data, args.output, args.source, args.prefixed)
