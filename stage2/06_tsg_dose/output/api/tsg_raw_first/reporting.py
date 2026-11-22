"""Render the new dose study report from its recorded summary CSV."""
import math

from stage2_util import palette, plotting

LABELS = {"pcc": "Pearson", "spearman": "Spearman", "rmse": "RMSE", "r2": "R²"}
SOURCE_COLORS = {"tenv": "blue", "tahoe_bulk": "orange"}


def render(summary_csv, output):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    rows = plotting.read_rows(summary_csv)
    if not rows:
        raise ValueError("raw-first summary CSV is empty")
    figure = Figure(figsize=(8, 5), layout="constrained")
    FigureCanvasAgg(figure)
    axes = figure.subplots(2, 2).flat
    labels = [f"{row['source']}\n{row['scenario']} · {row['fold']}" for row in rows]
    colors = [palette.color(SOURCE_COLORS[row["source"]], palette.FILL_SHADE) for row in rows]
    edges = [palette.color(SOURCE_COLORS[row["source"]], palette.LINE_SHADE) for row in rows]
    for axis, (metric, label) in zip(axes, LABELS.items(), strict=True):
        for index, row in enumerate(rows):
            value = float(row[metric]) if row[metric] else math.nan
            if math.isfinite(value):
                axis.bar(index, value, color=colors[index], edgecolor=edges[index])
            else:
                axis.scatter(index, 0, marker="x", color=palette.color("gray", palette.REFERENCE_SHADE))
        axis.set_xticks(range(len(rows)), labels, rotation=25, ha="right")
        axis.set_title(f"Δ {label}")
        axis.grid(axis="y", color=palette.color("gray", 8))
        axis.set_axisbelow(True)
    plotting.save_figure(figure, output, dpi=180)
