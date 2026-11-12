#!/usr/bin/env python3
"""Render the published gene-axis source-summary panel from its plotted-values CSV alone.

This is a publication script, not scientific code: it reads only
``figure/data/gene_axis_source_summary.csv`` and writes the PNG. It never opens
``output/`` products, matrices, checkpoints or truth, and it recomputes no
metric. The CSV carries every plotted value — the per-source HVG/DEG shares, the
fused-score histogram with the candidate-inflection marker (or the explicit
"no inflection" note), the tenv non-HVG curve coordinates and the source overlap
matrix cells — so this renderer needs no other input.

    python3 stage2/01_bases/figure/script/render_gene_axis_source_summary.py \
        --csv stage2/01_bases/figure/data/gene_axis_source_summary.csv \
        --png stage2/01_bases/figure/gene_axis_source_summary.png

Re-rendering from the CSV reproduces the committed PNG byte for byte.
"""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402  (backend must be set first)
import numpy as np  # noqa: E402

RECORDS = ("figure", "panel", "bar", "hist", "marker", "curve", "matrix")
PANELS = ("a", "b", "c", "d")
SHARES = ("hvg_share", "deg_up_share", "deg_down_share")
SHARE_LABELS = ("HVG share", "DEG up share", "DEG down share")
COLORS = ("#4C72B0", "#DD8452", "#55A868")


def load(path):
    """Read and validate the complete plotted-values table."""
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("plotted-values CSV is empty")
    kinds = {row["record"] for row in rows}
    if not set(RECORDS) <= kinds:
        raise ValueError(f"plotted-values CSV is incomplete: {sorted(kinds)}")
    figures = [row for row in rows if row["record"] == "figure"]
    if len(figures) != 1 or not figures[0]["title"]:
        raise ValueError("exactly one figure record with a title is required")
    panels = {row["panel"]: row["panel_title"] for row in rows if row["record"] == "panel"}
    if set(panels) != set(PANELS) or any(not value for value in panels.values()):
        raise ValueError("panel records must define panels a, b, c and d")
    bars = {}
    for row in rows:
        if row["record"] != "bar":
            continue
        if not row["source"]:
            raise ValueError("bar record without a source")
        bars[row["source"]] = {
            name: float(row[name]) if row[name] else None for name in SHARES}
    if not bars:
        raise ValueError("no per-source bar values")
    histogram = [(float(row["bin_left"]), float(row["bin_right"]), int(row["count"]))
                 for row in rows if row["record"] == "hist"]
    if not histogram:
        raise ValueError("fused-score histogram values are missing")
    markers = [row for row in rows if row["record"] == "marker"]
    if len(markers) != 1:
        raise ValueError("exactly one marker record is required")
    marker = {"n": int(markers[0]["marker_n"]) if markers[0]["marker_n"] else None,
              "score": float(markers[0]["marker_score"]) if markers[0]["marker_score"] else None,
              "note": markers[0]["note"]}
    curves = {}
    for row in rows:
        if row["record"] != "curve":
            continue
        key = int(row["k"])
        curves.setdefault(key, []).append((int(row["n"]), float(row["non_hvg_share"])))
    if not curves:
        raise ValueError("tenv non-HVG curve values are missing")
    matrix = {}
    order = []
    for row in rows:
        if row["record"] != "matrix":
            continue
        key = (row["row_source"], row["col_source"])
        matrix[key] = float(row["jaccard"])
        for source in key:
            if source not in order:
                order.append(source)
    if not matrix or len({row for pair in matrix for row in pair}) != len(order):
        raise ValueError("source overlap matrix is incomplete")
    if any((source, source) not in matrix for source in order):
        raise ValueError("source overlap matrix is missing a diagonal cell")
    return figures[0]["title"], panels, bars, histogram, marker, curves, matrix, order


def _short(name):
    """Human label for a CSV set name; no data access."""
    return (name.replace("_deg", "").replace("tenv_hvg_2000", "tenv HVG")
            .replace("dtoxs_ipsc", "dtoxs").replace("gse244740", "244740"))


def render(title, panels, bars, histogram, marker, curves, matrix, order, output):
    """Draw the four-panel summary with fixed publication geometry."""
    output = Path(output)
    if output.suffix != ".png":
        raise ValueError("a .png output path is required")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.6), layout="constrained")
    try:
        axis = axes[0][0]
        sources = list(bars)
        positions = np.arange(len(sources), dtype=float)
        width = 0.26
        for offset, (name, label, color) in enumerate(zip(SHARES, SHARE_LABELS, COLORS)):
            values = [bars[source][name] if bars[source][name] is not None else np.nan
                      for source in sources]
            axis.bar(positions + (offset - 1) * width, values, width=width,
                     color=color, label=label)
        axis.set_xticks(positions)
        axis.set_xticklabels([_short(source) for source in sources], rotation=35,
                             ha="right", fontsize=8)
        axis.set(ylabel="share", title=panels["a"])
        axis.grid(alpha=0.18, axis="y")
        axis.legend(frameon=False, fontsize=8)

        axis = axes[0][1]
        left = [bin_left for bin_left, _, _ in histogram]
        width = [bin_right - bin_left for bin_left, bin_right, _ in histogram]
        axis.bar(left, [count for _, _, count in histogram], width=width,
                 align="edge", color="#4C72B0")
        if marker["score"] is not None:
            axis.axvline(marker["score"], color="#C44E52", linewidth=1.8,
                         label=f"candidate N={marker['n']}")
            axis.legend(frameon=False, fontsize=8)
        axis.set(xlabel="fused rank_sum", ylabel="genes", title=panels["b"])
        axis.grid(alpha=0.18, axis="y")
        axis.text(0.02, 0.85, marker["note"], transform=axis.transAxes, va="top",
                  fontsize=7.5)

        axis = axes[1][0]
        for index, key in enumerate(sorted(curves)):
            points = sorted(curves[key])
            axis.plot([point[0] for point in points], [point[1] for point in points],
                      color=plt.get_cmap("viridis")(0.15 + 0.35 * index),
                      linewidth=1.6, label=f"tenv HVG K={key}")
        axis.set(xlabel="N (fused top-N)", ylabel="non-HVG share", title=panels["c"],
                 ylim=(0, 1.02))
        axis.grid(alpha=0.18)
        axis.legend(frameon=False, fontsize=8)

        axis = axes[1][1]
        grid = np.array([[matrix[(row, column)] for column in order] for row in order])
        image = axis.imshow(grid, cmap="viridis", vmin=0.0, vmax=1.0)
        axis.set_xticks(range(len(order)))
        axis.set_yticks(range(len(order)))
        axis.set_xticklabels([_short(source) for source in order], rotation=45,
                             ha="right", fontsize=7)
        axis.set_yticklabels([_short(source) for source in order], fontsize=7)
        for row in range(len(order)):
            for column in range(len(order)):
                value = grid[row, column]
                axis.text(column, row, f"{value:.2f}", ha="center", va="center",
                          fontsize=5.5,
                          color="white" if value < 0.55 else "black")
        axis.set_title(panels["d"])
        fig.colorbar(image, ax=axis, fraction=0.046, pad=0.02, label="Jaccard")

        fig.suptitle(title, fontsize=13)
        fig.savefig(output, dpi=160)
    finally:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", type=Path, required=True, help="figure/data plotted values")
    parser.add_argument("--png", type=Path, required=True, help="figure PNG to write")
    args = parser.parse_args()
    render(*load(args.csv), args.png)
    print(f"rendered {args.png} from {args.csv}")


if __name__ == "__main__":
    main()
