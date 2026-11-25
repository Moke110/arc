"""Render fold scores and mean paired differences from the Step 06.1 CSVs."""

import argparse
from itertools import product
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

STAGE = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(STAGE / "util/output/api"), str(STAGE / "05_tsg_loss/figure/script")]

from loss_plot_style import FOLD_HUES, metric_direction
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
FIELDS = ("value", "score_difference")
MEAN_HUE = 300
FIGURES = {"value": "delta-metrics.png", "score_difference": "delta-metrics-paired-difference.png"}


def read_table(path, keys):
    rows = {}
    for row in plotting.read_rows(path):
        key = tuple(int(row[name]) if name == "fold" else row[name] for name in keys)
        if key in rows:
            raise ValueError(f"duplicate plotted score in {path.name}: {key}")
        values = {field: float(row[field]) if row[field] else math.nan for field in FIELDS}
        if any(math.isinf(value) for value in values.values()):
            raise ValueError(f"infinite plotted score in {path.name}: {key}")
        rows[key] = values
    return rows


def read_scores(data):
    summary = read_table(data / "delta-metrics.csv", ("shape", "metric"))
    scores = read_table(data / "fold-scores.csv", ("shape", "fold", "metric"))
    shapes = tuple(dict.fromkeys(shape for shape, _ in summary))
    folds = tuple(sorted({fold for _, fold, _ in scores}))
    if (not shapes or shapes[0] != "widen_z-100" or folds not in ((0,), (0, 1, 2))
            or set(summary) != set(product(shapes, LABELS)) or set(scores) != set(product(shapes, folds, LABELS))):
        raise ValueError("shape plotting CSVs require a complete paired score grid")
    for (shape, fold, metric), score in scores.items():
        baseline = scores["widen_z-100", fold, metric]["value"]
        expected = score["value"] - baseline
        if not equal_score(score["score_difference"], expected):
            raise ValueError(f"paired difference disagrees with the control for {shape}, fold {fold}, {metric}")
    for (shape, metric), aggregate in summary.items():
        for field in FIELDS:
            values = [scores[shape, fold, metric][field] for fold in folds
                      if math.isfinite(scores[shape, fold, metric][field])]
            expected = math.fsum(values) / len(values) if values else math.nan
            if not equal_score(aggregate[field], expected):
                raise ValueError(f"mean disagrees with the plotted folds for {shape}, {metric}")
    return shapes, folds, scores, summary


def equal_score(left, right):
    return (math.isnan(left) and math.isnan(right)) or math.isclose(left, right, rel_tol=1e-13, abs_tol=1e-15)


def render(data, output):
    shapes, folds, scores, summary = read_scores(data)
    x = list(range(len(shapes)))
    offsets = {fold: (index - (len(folds) - 1) / 2) * 0.14 - 0.12 for index, fold in enumerate(folds)}
    for field, filename in FIGURES.items():
        difference = field == "score_difference"
        figure, plots = plt.subplots(3, 3, figsize=(15, 10), layout="constrained")
        for (metric, label), plot in zip(LABELS.items(), plots.flat, strict=True):
            metric_direction(plot, metric, label)
            for fold in folds:
                values = [scores[shape, fold, metric][field] for shape in shapes]
                plotting.scatter(plot, [position + offsets[fold] for position in x], values,
                                 hue=FOLD_HUES[fold], s=18, zorder=3)
            means = [summary[shape, metric][field] for shape in shapes]
            plotting.scatter(plot, [position + 0.21 for position in x], means,
                             hue=MEAN_HUE, marker="D", s=22, zorder=4)
            for position, mean in zip(x, means, strict=True):
                if math.isnan(mean):
                    plot.text(position, .02, "NA", ha="center", fontsize=7, transform=plot.get_xaxis_transform())
            if difference:
                plot.axhline(0, color="#666666", linewidth=0.8, zorder=2)
            plot.set_xticks(x, shapes, rotation=30, ha="right", fontsize=8)
            plot.set_xlim(-0.5, len(x) - 0.5)
            plot.grid(axis="y", color="#777777", alpha=0.25, zorder=1)
            plot.spines[["top", "right"]].set_visible(False)
        handles = [Line2D([], [], marker="o", linestyle="none", markersize=5,
                          color=palette.color_pair(FOLD_HUES[fold], 2).main, label=f"Fold {fold}") for fold in folds]
        handles.append(Line2D([], [], marker="D", linestyle="none", markersize=5,
                              color=palette.color_pair(MEAN_HUE, 2).main,
                              label="Equal-fold mean"))
        figure.legend(handles=handles, ncol=len(handles), loc="outside lower center", frameon=False)
        figure.supylabel("Score difference from the original architecture" if difference else "Held-out fold score")
        plotting.save_figure(figure, output / filename, dpi=180)
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == "__main__":
    main()
