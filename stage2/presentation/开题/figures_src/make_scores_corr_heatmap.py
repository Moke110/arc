#!/usr/bin/env python3
"""Slide-ready redraw of the cross-method Spearman rho heatmap (R11).

Reproduces the SAME 12x12 matrix as ``plot_scores_spearman`` in
``stage1/xpert-fda-miri-v2/09_miri_efficacy_scores/miri_efficacy_scores.py``
from that step's ``output/scores.csv`` (N=1004 drugs), but with LARGER fonts and
the long 3-line title removed.

    python3 presentation/figures_src/make_scores_corr_heatmap.py

Output: presentation/figures/scores_corr_heatmap.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr

ROOT = Path("/home/chang/arc")
SCORES = ROOT / "stage1/xpert-fda-miri-v2/09_miri_efficacy_scores/output/scores.csv"
OUT = ROOT / "presentation/figures/scores_corr_heatmap.png"

CRITICAL = [
    "GO-Inflammatory Response",
    "GO-Regulation Of Inflammatory Response",
    "KEGG-Apoptosis",
    "GO-Cellular Response To Oxidative Stress",
    "KEGG-Necroptosis",
    "GO-Regulation Of Cell Death",
    "KEGG-Calcium signaling pathway",
    "GO-Mitochondrial Gene Expression",
]

SCORE_COLS = ["WCS_zs", "RGES_zs", "Zhang_zs", "pathway_reversal_zs"]
LABELS = [
    "\u2212WCS", "\u2212RGES", "\u2212Zhang", "Pathway Reversal",
    "Inflammation", "Inflamm. Reg.", "Apoptosis", "Oxid. Stress",
    "Necroptosis", "Cell Death Reg.", "Ca Signaling", "Mito. Gene",
]
NEGATE = {"WCS_zs", "RGES_zs", "Zhang_zs"}


def main():
    scores = pd.read_csv(SCORES)
    n_drugs = len(scores)
    pw_cols = [f"{n}_efficacy" for n in CRITICAL]
    all_cols = SCORE_COLS + pw_cols

    n = len(all_cols)
    data_mat = np.zeros((n_drugs, n))
    for j, col in enumerate(all_cols):
        vals = scores[col].values.astype(float)
        if col in NEGATE:
            vals = -vals
        data_mat[:, j] = vals

    rho = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            rho[i, j], _ = spearmanr(data_mat[:, i], data_mat[:, j])

    w2b = LinearSegmentedColormap.from_list("w2b", ["#ffffff", "#2166ac"], N=256)
    fig, ax = plt.subplots(figsize=(11, 10))
    im = ax.imshow(rho, cmap=w2b, vmin=0, vmax=1, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_xticklabels(LABELS, fontsize=15, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(LABELS, fontsize=15)
    for i in range(n):
        for j in range(n):
            color = "white" if rho[i, j] > 0.5 else "black"
            ax.text(j, i, f"{rho[i, j]:.2f}", ha="center", va="center",
                    fontsize=12, fontweight="bold", color=color)
    cbar = plt.colorbar(im, ax=ax, shrink=0.82)
    cbar.set_label("Spearman \u03c1", fontsize=15)
    cbar.ax.tick_params(labelsize=14)
    ax.set_title("Spearman \u03c1", fontsize=16, fontweight="bold", pad=12)
    fig.tight_layout()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT}  (N={n_drugs})")

    def v(a, b):
        return rho[LABELS.index(a), LABELS.index(b)]

    print(f"diagonal = {rho[0, 0]:.2f}")
    print(f"\u2212Zhang vs Apoptosis = {v('\u2212Zhang', 'Apoptosis'):.2f}")
    print(f"\u2212WCS vs \u2212Zhang = {v('\u2212WCS', '\u2212Zhang'):.2f}")
    print(f"\u2212Zhang vs Oxid. Stress = {v('\u2212Zhang', 'Oxid. Stress'):.2f}")
    print(f"Ca Signaling vs Mito. Gene = {v('Ca Signaling', 'Mito. Gene'):.2f}")


if __name__ == "__main__":
    main()
