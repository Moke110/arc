#!/usr/bin/env python3
"""Redraw the MIRI filter-progression funnel as a presentation figure (R10).

Presentation copy of `10_filter/filter.py::plot_pyramid` for 开题 slide 20:
larger fonts, no background grid, and no title. Counts are read from the real
workflow outputs (row counts), never hard-coded.

Output (PNG, 200 dpi) -> presentation/figures/miri_filter_progression.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path("/home/chang/arc")
SCORES = ROOT / "stage1/xpert-fda-miri-v2/09_miri_efficacy_scores/output/scores.csv"
WCS = ROOT / "stage1/xpert-fda-miri-v2/10_filter/output/wcs_filtered.csv"
PATHF = ROOT / "stage1/xpert-fda-miri-v2/10_filter/output/path_filtered.csv"
OUT = ROOT / "presentation/figures/miri_filter_progression.png"

plt.rcParams.update({
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "font.size": 22,
    "axes.labelsize": 24,
    "ytick.labelsize": 22,
    "font.family": "DejaVu Sans",
})


def n_rows(path):
    return len(pd.read_csv(path))


def main():
    n_all = n_rows(SCORES)
    n_wcs = n_rows(WCS)
    n_path = n_rows(PATHF)
    print(f"All FDA Drugs: {n_all}")
    print(f"WCS < 0:      {n_wcs}")
    print(f"Pathway >=6/8: {n_path}")

    stages = [
        ("All FDA Drugs", n_all, "#607d8b"),
        ("WCS < 0", n_wcs, "#2196f3"),
        ("Pathway >= 6/8", n_path, "#4caf50"),
    ]
    counts = [s[1] for s in stages]

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.barh(range(len(stages)), counts, height=0.6,
                   color=[s[2] for s in stages], edgecolor="white", linewidth=2, alpha=0.9)
    max_w = max(counts)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_width() + max_w * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{count} ({count / n_all * 100:.1f}%)",
                va="center", fontsize=20, fontweight="bold", color="#333333")
    ax.set_yticks(range(len(stages)))
    ax.set_yticklabels([s[0] for s in stages], fontsize=22)
    ax.set_xlim(0, max_w * 1.2)
    ax.invert_yaxis()
    ax.set_xlabel("Number of Drugs", fontsize=24)
    ax.grid(False)
    ax.set_title("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Saved:", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
