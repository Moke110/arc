#!/usr/bin/env python3
"""Redraw the MIRI efficacy-score distributions as standalone figures.

For the 开题 slides 15-17 【补】 tasks: each score distribution is redrawn on its
own, with the background grid and the figure title removed and the fonts
enlarged. Data: stage1/xpert-fda-miri-v3/09_miri_efficacy_scores/output/scores.csv

Outputs (PNG, 200 dpi) -> presentation/figures/
    wcs_distribution.png
    rges_distribution.png
    zhang_distribution.png
    pathway_reversal_distribution.png
    pathway_<short>.png   (the 8 critical-pathway efficacies)
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/home/chang/arc")
# The 开题 deck's efficacy-score figures (hist_scores.png) were produced by the
# v2 workflow; redraw from the same table so the numbers stay consistent.
SCORES = ROOT / "stage1/xpert-fda-miri-v2/09_miri_efficacy_scores/output/scores.csv"
OUT = ROOT / "presentation/figures"
OUT.mkdir(parents=True, exist_ok=True)

N_PATHS = 40  # total MIRI pathways (pathway reversal denominator)

CRITICAL = [
    ("GO-Inflammatory Response", "Inflammation", "firebrick"),
    ("GO-Regulation Of Inflammatory Response", "Inflamm. Reg.", "firebrick"),
    ("KEGG-Apoptosis", "Apoptosis", "firebrick"),
    ("GO-Cellular Response To Oxidative Stress", "Oxid. Stress", "firebrick"),
    ("KEGG-Necroptosis", "Necroptosis", "firebrick"),
    ("GO-Regulation Of Cell Death", "Cell Death Reg.", "firebrick"),
    ("KEGG-Calcium signaling pathway", "Ca Signaling", "steelblue"),
    ("GO-Mitochondrial Gene Expression", "Mito. Gene", "steelblue"),
]

plt.rcParams.update({
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "font.size": 17,
    "axes.labelsize": 19,
    "axes.titlesize": 19,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 14,
    "axes.linewidth": 1.5,
    "font.family": "DejaVu Sans",
})


def clean(ax):
    ax.grid(False)
    ax.set_title("")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def save(fig, name):
    out = OUT / name
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Saved:", out.relative_to(ROOT))


def hist(ax, vals, bins, color, xlabel, pct=False):
    weights = None
    if pct:
        weights = np.ones(len(vals)) * 100.0 / len(vals)
    ax.hist(vals, bins=bins, color=color, alpha=0.8, edgecolor="white",
            linewidth=0.4, weights=weights)
    ax.axvline(0, color="black", linestyle="--", linewidth=1.2, alpha=0.6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Percentage (%)" if pct else "Count")


def main():
    df = pd.read_csv(SCORES)
    print(f"Loaded {len(df)} drugs")

    # --- WCS ---
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    v = df["WCS"].values
    hist(ax, v, np.linspace(v.min(), v.max(), 60), "steelblue", "WCS")
    clean(ax)
    save(fig, "wcs_distribution.png")

    # --- RGES ---
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    hist(ax, df["RGES"].values, 60, "seagreen", "RGES")
    clean(ax)
    save(fig, "rges_distribution.png")

    # --- Zhang ---
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    hist(ax, df["Zhang_c"].values, 60, "coral", "Zhang c")
    clean(ax)
    save(fig, "zhang_distribution.png")

    # --- Pathway reversal (fraction of 40 pathways reversed) ---
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    v = df["pathway_reversal"].values
    step = 1.0 / N_PATHS
    bins = np.arange(v.min() - step / 2, v.max() + step, step)
    ax.hist(v, bins=bins, color="mediumpurple", alpha=0.8, edgecolor="white",
            linewidth=0.5)
    ax.set_xlabel("Pathway reversal rate")
    ax.set_ylabel("Count")
    clean(ax)
    save(fig, "pathway_reversal_distribution.png")

    # --- 8 critical pathways ---
    bins = np.arange(-1, 1 + 0.02, 0.02)
    for name, short, color in CRITICAL:
        vals = df[f"{name}_efficacy"].values
        fig, ax = plt.subplots(figsize=(5.6, 4.3))
        hist(ax, vals, bins, color, short)
        clean(ax)
        slug = short.replace(".", "").replace(" ", "_")
        save(fig, f"pathway_{slug}.png")


if __name__ == "__main__":
    main()
