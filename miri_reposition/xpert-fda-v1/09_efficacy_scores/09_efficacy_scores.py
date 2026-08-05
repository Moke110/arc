#!/usr/bin/env python3
"""
Step 09 — Efficacy scores.

Runs the three MIRI efficacy scorers (WCS+RGES, Zhang/ssCMap, pathway RGES)
against the integrated MIRI reference (08_miri_ref_path/output/) and the
curated pathway gene sets (input/MIRI_pathways.json), producing the per-drug
z-component table for all 3,127 drugs. The integrated score is NOT computed
here — step 11 integrates after filtering.

Inputs:
  07_zscore/output/xpert_fda-cm_deg-zs.csv   (z-scored predicted DEG)
  08_miri_ref_path/output/MIRI_up-down.json
  08_miri_ref_path/output/MIRI_rank.json
  input/MIRI_pathways.json

Outputs:
  09_efficacy_scores/output/scores.csv       — raw scores + z-components
  09_efficacy_scores/hist_scores.png             — raw score distributions (2x3)
  09_efficacy_scores/hist_critical-path-efficacy.png — per-critical-pathway efficacies
  09_efficacy_scores/heatmap_scores-spearman.png — cross-method Spearman rho
  09_efficacy_scores/hist_stratified.png         — WCS + 8 critical pathways x 3 methods

Usage:
    python 09_efficacy_scores/09_efficacy_scores.py
    python 09_efficacy_scores/09_efficacy_scores.py --skip-scores   # replot only
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT  # noqa: E402
from scoring import (  # noqa: E402
    load_critical_names,
    merge_components,
    score_all_pathways,
    score_all_wcs_rges,
    score_all_zhang,
)

STEP_DIR = WORKFLOW_ROOT / "09_efficacy_scores"
DEG = WORKFLOW_ROOT / "07_zscore" / "output" / "xpert_fda-cm_deg-zs.csv"
GENE_SET = WORKFLOW_ROOT / "08_miri_ref_path" / "output" / "MIRI_up-down.json"
RANK_REF = WORKFLOW_ROOT / "08_miri_ref_path" / "output" / "MIRI_rank.json"
PATHWAYS = WORKFLOW_ROOT / "input" / "MIRI_pathways.json"
OUTPUT = STEP_DIR / "output" / "scores.csv"

CRITICAL_SHORT = {
    "GO-Inflammatory Response": "Inflammation",
    "GO-Regulation Of Inflammatory Response": "Inflamm. Reg.",
    "KEGG-Apoptosis": "Apoptosis",
    "GO-Cellular Response To Oxidative Stress": "Oxid. Stress",
    "KEGG-Necroptosis": "Necroptosis",
    "GO-Regulation Of Cell Death": "Cell Death Reg.",
    "KEGG-Calcium signaling pathway": "Ca Signaling",
    "GO-Mitochondrial Gene Expression": "Mito. Gene",
}


def load_pathway_meta():
    with open(PATHWAYS) as f:
        meta = json.load(f)
    critical = [p["name"] for p in meta["pathways"] if p.get("critical")]
    return meta, critical


def save(fig, name):
    out = STEP_DIR / name
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.relative_to(WORKFLOW_ROOT)}")


def plot_critical_path_efficacy(scores, pw_meta, critical):
    bins = np.arange(-1, 1 + 0.02, 0.02)
    fig, axes = plt.subplots(4, 2, figsize=(14, 20))
    axes_flat = axes.flatten()
    for idx, name in enumerate(critical):
        ax = axes_flat[idx]
        p = next(p for p in pw_meta["pathways"] if p["name"] == name)
        vals = scores[f"{name}_efficacy"].values
        vals_mean = vals.mean()
        color = "firebrick" if p["direction"] == "UP" else "steelblue"
        ax.hist(vals, bins=bins, color=color, alpha=0.75, edgecolor="white", linewidth=0.2)
        ax.axvline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.axvline(vals_mean, color="darkred" if p["direction"] == "UP" else "darkblue",
                   linestyle="-", linewidth=1.2, label=f"Mean = {vals_mean:+.4f}")
        ax.set_xlabel("Efficacy (signed RGES)")
        ax.set_ylabel("Count")
        ax.set_xlim(-1, 1)
        ax.set_title(f"{CRITICAL_SHORT[name]}  ({p['direction']}, {p['n_genes']} genes)",
                     fontsize=11)
        ax.legend(fontsize=8)
    for idx in range(len(critical), 8):
        axes_flat[idx].set_visible(False)
    fig.tight_layout()
    save(fig, "hist_critical-path-efficacy.png")


def plot_hist_scores(scores, pw_meta, critical):
    n_paths = len(pw_meta["pathways"])
    crit_up = [n for n in critical
               if any(p["name"] == n and p["direction"] == "UP" for p in pw_meta["pathways"])]
    crit_down = [n for n in critical
                 if any(p["name"] == n and p["direction"] == "DOWN" for p in pw_meta["pathways"])]
    all_up_eff = scores[[f"{n}_efficacy" for n in crit_up]].values.ravel()
    all_down_eff = scores[[f"{n}_efficacy" for n in crit_down]].values.ravel()

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    ax = axes[0, 0]
    wcs_vals = scores["WCS"].values
    nonzero_wcs = wcs_vals[wcs_vals != 0]
    ax.hist(wcs_vals, bins=np.linspace(wcs_vals.min(), wcs_vals.max(), 60),
            color="steelblue", alpha=0.7)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("WCS")
    ax.set_ylabel("Count")
    ax.set_title(f"WCS  |  mean={wcs_vals.mean():+.4f}  std={wcs_vals.std():.4f}  |  "
                 f"Nonzero: {len(nonzero_wcs)}/{len(wcs_vals)}")

    ax = axes[0, 1]
    ax.hist(scores["RGES"], bins=60, color="seagreen", alpha=0.7)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("RGES")
    ax.set_ylabel("Count")
    ax.set_title(f"RGES  |  mean={scores['RGES'].mean():+.4f}  std={scores['RGES'].std():.4f}")

    ax = axes[0, 2]
    ax.hist(all_up_eff, bins=60, color="firebrick", alpha=0.7,
            weights=np.ones(len(all_up_eff)) * 100 / len(all_up_eff))
    ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(all_up_eff.mean(), color="darkred", linestyle="-", linewidth=1.2,
               label=f"Mean = {all_up_eff.mean():+.4f}")
    ax.set_xlabel("MIRI-UP Pathway Efficacy")
    ax.set_ylabel("Percentage (%)")
    ax.set_title(f"UP Pathways ({len(crit_up)})  |  mean={all_up_eff.mean():+.4f}  "
                 f"std={all_up_eff.std():.4f}  N={len(all_up_eff)}")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.hist(scores["Zhang_c"], bins=60, color="coral", alpha=0.7)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Zhang c")
    ax.set_ylabel("Count")
    strong = (scores["Zhang_c"] < -0.3).sum()
    ax.set_title(f"Zhang c  |  mean={scores['Zhang_c'].mean():+.4f}  "
                 f"std={scores['Zhang_c'].std():.4f}  |  c<-0.3: {strong}")

    ax = axes[1, 1]
    pr_vals = scores["pathway_reversal"].values
    pr_step = 1.0 / n_paths
    pr_bins = np.arange(pr_vals.min() - pr_step / 2, pr_vals.max() + pr_step, pr_step)
    ax.hist(pr_vals, bins=pr_bins, color="mediumpurple", alpha=0.7,
            edgecolor="white", linewidth=0.5)
    ax.axvline(pr_vals.mean(), color="indigo", linestyle="--", linewidth=0.8,
               label=f"Mean = {pr_vals.mean():.4f}")
    ax.set_xlabel("Pathway Reversal")
    ax.set_ylabel("Count")
    ax.set_title(f"Pathway Reversal  |  mean={pr_vals.mean():.4f}  std={pr_vals.std():.4f}")
    ax.legend(fontsize=9)

    ax = axes[1, 2]
    ax.hist(all_down_eff, bins=60, color="steelblue", alpha=0.7,
            weights=np.ones(len(all_down_eff)) * 100 / len(all_down_eff))
    ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(all_down_eff.mean(), color="darkblue", linestyle="-", linewidth=1.2,
               label=f"Mean = {all_down_eff.mean():+.4f}")
    ax.set_xlabel("MIRI-DOWN Pathway Efficacy")
    ax.set_ylabel("Percentage (%)")
    ax.set_title(f"DOWN Pathways ({len(crit_down)})  |  mean={all_down_eff.mean():+.4f}  "
                 f"std={all_down_eff.std():.4f}  N={len(all_down_eff)}")
    ax.legend(fontsize=8)

    fig.tight_layout()
    save(fig, "hist_scores.png")


def plot_scores_spearman(scores, critical):
    n_drugs = len(scores)
    score_cols = ["WCS_zs", "RGES_zs", "Zhang_zs", "pathway_reversal_zs"]
    score_labels = ["−WCS", "−RGES", "−Zhang", "Pathway\nReversal"]
    negate_cols = {"WCS_zs", "RGES_zs", "Zhang_zs"}
    pw_cols = [f"{n}_efficacy" for n in critical]
    pw_labels = [CRITICAL_SHORT[n] for n in critical]
    all_cols = score_cols + pw_cols
    all_labels = score_labels + pw_labels

    n = len(all_cols)
    data_mat = np.zeros((n_drugs, n))
    for j, col in enumerate(all_cols):
        vals = scores[col].values.astype(float)
        if col in negate_cols:
            vals = -vals
        data_mat[:, j] = vals

    rho_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            rho_matrix[i, j], _ = spearmanr(data_mat[:, i], data_mat[:, j])

    w2b = LinearSegmentedColormap.from_list("w2b", ["#ffffff", "#2166ac"], N=256)
    fig, ax = plt.subplots(figsize=(12, 10.5))
    im = ax.imshow(rho_matrix, cmap=w2b, vmin=0, vmax=1, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_xticklabels(all_labels, fontsize=7, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(all_labels, fontsize=7)
    for i in range(n):
        for j in range(n):
            color = "white" if rho_matrix[i, j] > 0.5 else "black"
            ax.text(j, i, f"{rho_matrix[i, j]:.2f}", ha="center", va="center",
                    fontsize=6.5, fontweight="bold", color=color)
    cbar = plt.colorbar(im, ax=ax, shrink=0.82)
    cbar.set_label("Spearman ρ", fontsize=12)
    ax.set_title(f"Cross-Method Spearman ρ (N={n_drugs})\n"
                 "4 Z-scores (direction-aligned) + 8 Critical Pathway Efficacy\n"
                 "−WCS, −RGES, −Zhang: negated so higher = better therapeutic effect",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    save(fig, "heatmap_scores-spearman.png")


def plot_hist_stratified(scores, pw_meta, critical):
    n_paths = len(pw_meta["pathways"])
    method_cols = ["RGES", "Zhang_c", "pathway_reversal"]
    method_labels = ["RGES", "Zhang", "Pathway Reversal"]
    row_labels = ["WCS"] + [CRITICAL_SHORT[n] for n in critical]

    fig, axes = plt.subplots(9, 3, figsize=(18, 45))
    for row_idx, row_label in enumerate(row_labels):
        if row_idx == 0:
            wcs_vals = scores["WCS"].values
            masks = [(wcs_vals < 0, "#2166ac", "WCS < 0"),
                     (wcs_vals == 0, "#d9d9d9", "WCS = 0"),
                     (wcs_vals > 0, "#b2182b", "WCS > 0")]
            for col_idx, (label, col) in enumerate(zip(method_labels, method_cols)):
                ax = axes[row_idx, col_idx]
                vals = scores[col].values
                if col_idx == 2:
                    pr_step = 1.0 / n_paths
                    bins = np.arange(-0.5 * pr_step, 1.0 + 1.5 * pr_step, pr_step)
                else:
                    bins = np.histogram_bin_edges(vals, bins=50)
                bc = (bins[:-1] + bins[1:]) / 2
                bw = bins[1] - bins[0]
                bottom = np.zeros(len(bc))
                for mask, color, mlabel in masks:
                    h, _ = np.histogram(vals[mask], bins=bins)
                    ax.bar(bc, h, width=bw, bottom=bottom, color=color, alpha=0.9,
                           edgecolor="none", label=f"{mlabel} ({mask.sum()})")
                    bottom += h
                ax.axvline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
                ax.set_xlabel(label)
                if col_idx == 0:
                    ax.set_ylabel("WCS\nCount")
                ax.set_title(f"{label} by WCS Category", fontsize=11, fontweight="bold")
                ax.legend(fontsize=7)
        else:
            pw_name = critical[row_idx - 1]
            p_meta = next(p for p in pw_meta["pathways"] if p["name"] == pw_name)
            direction = p_meta["direction"]
            eff_vals = scores[f"{pw_name}_efficacy"].values
            dark_color, light_color = "#2166ac", "#d1e5f0"
            mask_pos = eff_vals > 0
            mask_neg = eff_vals <= 0
            for col_idx, (label, col) in enumerate(zip(method_labels, method_cols)):
                ax = axes[row_idx, col_idx]
                vals = scores[col].values
                if col_idx == 2:
                    pr_step = 1.0 / n_paths
                    bins = np.arange(-0.5 * pr_step, 1.0 + 1.5 * pr_step, pr_step)
                else:
                    bins = np.histogram_bin_edges(vals, bins=50)
                h_pos, _ = np.histogram(vals[mask_pos], bins=bins)
                h_neg, _ = np.histogram(vals[mask_neg], bins=bins)
                bc = (bins[:-1] + bins[1:]) / 2
                bw = bins[1] - bins[0]
                ax.bar(bc, h_neg, width=bw, color=light_color, alpha=0.95,
                       edgecolor="none", label=f"Eff ≤ 0 ({mask_neg.sum()})")
                ax.bar(bc, h_pos, width=bw, bottom=h_neg, color=dark_color, alpha=0.95,
                       edgecolor="none", label=f"Eff > 0 ({mask_pos.sum()})")
                ax.axvline(0, color="black", linestyle="--", linewidth=0.7, alpha=0.4)
                ax.set_xlabel(label)
                if col_idx == 0:
                    ax.set_ylabel(f"{CRITICAL_SHORT[pw_name]}\n({direction})\nCount")
                ax.set_title(f"{label} by {CRITICAL_SHORT[pw_name]} Efficacy",
                             fontsize=11, fontweight="bold")
                ax.legend(fontsize=7)

    fig.suptitle(f"Stratified Score Distributions — WCS + 8 Critical Pathways (N={len(scores)})",
                 fontsize=15, fontweight="bold", y=1.005)
    fig.tight_layout()
    save(fig, "hist_stratified.png")


def main():
    parser = argparse.ArgumentParser(description="Compute MIRI efficacy score components")
    parser.add_argument("--deg-csv", default=str(DEG))
    parser.add_argument("--gene-set", default=str(GENE_SET))
    parser.add_argument("--rank-ref", default=str(RANK_REF))
    parser.add_argument("--pathways", default=str(PATHWAYS))
    parser.add_argument("--output", default=str(OUTPUT))
    parser.add_argument("--skip-scores", action="store_true",
                        help="Skip scoring; only regenerate plots from the existing table")
    args = parser.parse_args()

    sns.set_theme(style="whitegrid", context="notebook", palette="muted")
    plt.rcParams.update({"figure.dpi": 120, "figure.figsize": (10, 6),
                         "font.size": 11, "axes.titlesize": 14, "axes.labelsize": 12})

    pw_meta, critical = load_pathway_meta()

    if args.skip_scores:
        scores = pd.read_csv(args.output)
        print(f"Loaded cached scores: {len(scores)} drugs")
    else:
        # Input DEG is already z-scored (step 07 output) — do not re-normalize.
        wcs = score_all_wcs_rges(args.deg_csv, args.gene_set, zscore=False)
        zhang = score_all_zhang(args.deg_csv, args.rank_ref, zscore=False)
        pw, critical = score_all_pathways(args.deg_csv, args.pathways, zscore=False)
        scores = merge_components(wcs, zhang, pw, critical)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        scores.to_csv(args.output, index=False)
        print(f"Scores: {len(scores)} drugs x {len(scores.columns)} columns -> {args.output}")

    print("Score-distribution plots:")
    plot_critical_path_efficacy(scores, pw_meta, critical)
    plot_hist_scores(scores, pw_meta, critical)
    plot_scores_spearman(scores, critical)
    plot_hist_stratified(scores, pw_meta, critical)

    print(f"\nStep 09 complete — {STEP_DIR.relative_to(WORKFLOW_ROOT)}/")


if __name__ == "__main__":
    main()
