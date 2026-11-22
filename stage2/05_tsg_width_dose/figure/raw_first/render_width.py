"""Render raw-first width/dose test coordinates from summary CSV."""

import math

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from stage2_util import palette, plotting

LABELS = {"pcc": "Pearson r", "spearman": "Spearman ρ", "rmse": "RMSE", "r2": "R²"}
METHOD_COLORS = {
    "no_dose": "blue", "latent_numeric": "orange", "latent_embedding": "green", "FiLM": "purple",
}
WIDTHS = ("original", "1.5x", "2x", "2.5x")
FOLD_MARKERS = {"0": "o", "1": "s", "2": "^"}


def render(summary, output):
    rows = plotting.read_rows(summary)
    coordinates = [(row["width"], row["method"], row["fold"]) for row in rows]
    if not coordinates or len(set(coordinates)) != len(coordinates):
        raise ValueError("expected unique width, encoding and fold coordinates")
    if len({row["input_manifest_sha256"] for row in rows}) != 1:
        raise ValueError("summary mixes raw-first input generations")

    figure = Figure(figsize=(10, 5), layout="constrained")
    FigureCanvasAgg(figure)
    methods = tuple(METHOD_COLORS)
    axes = figure.subplots(2, 2).flat
    colors = {method: palette.color(color, palette.LINE_SHADE)
              for method, color in METHOD_COLORS.items()}
    for axis, (metric, title) in zip(axes, LABELS.items(), strict=True):
        for row, (width, method, fold) in zip(rows, coordinates, strict=True):
            value = float(row[metric]) if row[metric] else math.nan
            x = WIDTHS.index(width) + 0.18 * (methods.index(method) - 1.5)
            if math.isfinite(value):
                axis.scatter(x, value, marker=FOLD_MARKERS[fold], color=colors[method])
        axis.set_xlim(-0.5, len(WIDTHS) - 0.5)
        axis.set_xticks(range(len(WIDTHS)), WIDTHS)
        axis.set_xlabel("TranSiGen width")
        axis.set_title(f"Δ {title}")
        axis.grid(axis="y", color=palette.color("gray", 8))
        axis.set_axisbelow(True)

    legend = figure.axes[0]
    for method in methods:
        legend.scatter([], [], color=colors[method], label=method.replace("_", " "))
    for fold, marker in FOLD_MARKERS.items():
        legend.scatter([], [], color=palette.color("gray", palette.REFERENCE_SHADE),
                       marker=marker, label=f"Fold {fold}")
    handles, labels = legend.get_legend_handles_labels()
    figure.legend(handles, labels, loc="outside upper center", ncols=len(labels), frameon=False)
    plotting.save_figure(figure, output, dpi=180)
