#!/usr/bin/env python3
"""
Step 11 — Integrate.

Computes the integrated score over ALL drugs (z-components are population
statistics, so integration needs the full 3,127-drug table), then ranks the
filtered candidates from step 10 and writes the final report.

  integrated = -0.3*RGES_zs - 0.2*Zhang_zs + 0.3*pathway_reversal_zs
               + 0.025*sum(critical efficacy zs)

Inputs:
  09_efficacy_scores/output/scores.csv
  10_filter/output/path_filtered.csv
  input/MIRI_pathways.json

Outputs:
  11_integrate/output/integrated_scores.csv — all drugs + integrated score, sorted
  11_integrate/output/ranked_drugs.csv      — filtered drugs ranked by integrated
  11_integrate/output/report.md             — deterministic final report
  11_integrate/heatmap_rank.png             — score-percentile heatmap (all filtered)
  11_integrate/heatmap_rank@50.png          — same, integrated pct >= 0.5

Usage:
    python 11_integrate/11_integrate.py
"""

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import percentileofscore

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT  # noqa: E402
from scoring import add_integrated, load_critical_names  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "11_integrate"
SCORES = WORKFLOW_ROOT / "09_efficacy_scores" / "output" / "scores.csv"
PATH_FILTERED = WORKFLOW_ROOT / "10_filter" / "output" / "path_filtered.csv"
PATHWAYS = WORKFLOW_ROOT / "input" / "MIRI_pathways.json"

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

RWG = LinearSegmentedColormap.from_list("rwg", ["#d73027", "#f5f5f5", "#1a9850"], N=256)
LOWER_BETTER = {"WCS", "RGES", "Zhang_c"}


def therapeutic_percentile_matrix(frame, scores, heat_cols):
    """Percentile of each drug's value vs all drugs; 1.0 = therapeutic."""
    mat = np.zeros((len(frame), len(heat_cols)))
    for j, col in enumerate(heat_cols):
        full_vals = scores[col].values
        lb = col in LOWER_BETTER
        for i in range(len(frame)):
            pct = percentileofscore(full_vals, frame[col].values[i], kind="rank") / 100.0
            if not lb:
                pct = 1.0 - pct
            mat[i, j] = pct
    return 1.0 - mat


def render_rank_heatmap(merged, scores, critical, title, out_name):
    score_heat_cols = (["integrated", "WCS", "RGES", "Zhang_c", "pathway_reversal"]
                       + [f"{n}_efficacy" for n in critical])
    score_labels = (["Integrated", "WCS", "RGES", "Zhang", "Pathway Rev"]
                    + [CRITICAL_SHORT[n] for n in critical])

    merged = merged.sort_values("integrated", ascending=False).reset_index(drop=True)
    score_norm = therapeutic_percentile_matrix(merged, scores, score_heat_cols)

    drug_labels = merged["drug_name"].values
    n_rows = len(merged)
    n_score_cols = len(score_heat_cols)

    row_h, cell_w_score, cbar_w = 0.30, 0.42, 0.25
    hm_h = n_rows * row_h
    hm_w_score = n_score_cols * cell_w_score
    margin_left, margin_right = 3.2, 0.5
    margin_top, margin_bottom = 0.8, 0.7

    fig_w = margin_left + hm_w_score + 0.3 + cbar_w + margin_right
    fig_h = margin_top + hm_h + margin_bottom
    fig = plt.figure(figsize=(fig_w, fig_h))

    ax_score = fig.add_axes([margin_left / fig_w, margin_bottom / fig_h,
                             hm_w_score / fig_w, hm_h / fig_h])
    im_score = ax_score.imshow(score_norm, aspect="auto", cmap=RWG, vmin=0, vmax=1)
    ax_score.set_xticks(range(n_score_cols))
    ax_score.set_xticklabels(score_labels, rotation=45, ha="right", fontsize=12)
    ax_score.set_yticks(range(n_rows))
    ax_score.set_yticklabels(drug_labels, fontsize=12)
    ax_score.set_title("Therapeutic Percentile", fontsize=12, fontweight="bold", pad=6)

    cax_left = margin_left + hm_w_score + 0.3
    cax = fig.add_axes([cax_left / fig_w, (margin_bottom + hm_h * 0.55) / fig_h,
                        cbar_w / fig_w, (hm_h * 0.45) / fig_h])
    cb = plt.colorbar(im_score, cax=cax)
    cb.ax.tick_params(labelsize=12, pad=1)
    cb.set_label("Percentile in 3127 drugs", fontsize=12, fontweight="bold")
    cb.set_ticks([0, 0.5, 1])
    cb.set_ticklabels(["0\n(harmful)", "0.5", "1\n(therapeutic)"])

    fig.text(0.5, (fig_h - 0.10) / fig_h, title,
             fontsize=18, fontweight="bold", ha="center", va="top")

    out = STEP_DIR / out_name
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.relative_to(WORKFLOW_ROOT)}")


def write_report(ranked, n_all, n_wcs, integrated_full):
    pcts = [percentileofscore(integrated_full, v, kind="rank") / 100.0
            for v in ranked["integrated"].values]
    lines = [
        "# MIRI Drug Repositioning — Final Report",
        "",
        f"Generated: {date.today().isoformat()} (deterministic, step 11)",
        "",
        "## Screening funnel",
        "",
        f"- FDA-approved drugs screened: **{n_all}**",
        f"- WCS < 0: **{n_wcs}**",
        f"- WCS < 0 + >= 6/8 critical pathways reversed: **{len(ranked)}**",
        "",
        "## Integrated score",
        "",
        "`integrated = -0.3*RGES_zs - 0.2*Zhang_zs + 0.3*pathway_reversal_zs "
        "+ 0.025*sum(critical efficacy zs)`",
        "",
        "Computed over all drugs; candidates ranked by integrated score descending.",
        "",
        f"## Ranked candidates ({len(ranked)})",
        "",
        "| Rank | Drug | Integrated | Percentile |",
        "|---:|---|---:|---:|",
    ]
    for i, ((_, row), pct) in enumerate(zip(ranked.iterrows(), pcts), start=1):
        lines.append(f"| {i} | {row['drug_name']} | {row['integrated']:+.4f} | {pct:.3f} |")
    out = STEP_DIR / "output" / "report.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"  Saved: {out.relative_to(WORKFLOW_ROOT)}")


def main():
    parser = argparse.ArgumentParser(description="Integrate scores, rank candidates, report")
    parser.add_argument("--scores", default=str(SCORES))
    parser.add_argument("--path-filtered", default=str(PATH_FILTERED))
    parser.add_argument("--pathways", default=str(PATHWAYS))
    args = parser.parse_args()

    sns.set_theme(style="whitegrid", context="notebook", palette="muted")

    critical = load_critical_names(args.pathways)

    scores = pd.read_csv(args.scores)
    integrated = add_integrated(scores, critical).sort_values(
        "integrated", ascending=False).reset_index(drop=True)

    out_dir = STEP_DIR / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    integrated.to_csv(out_dir / "integrated_scores.csv", index=False)
    print(f"Integrated scores: {len(integrated)} drugs -> {out_dir / 'integrated_scores.csv'}")

    path_filtered = pd.read_csv(args.path_filtered)
    ranked = (path_filtered[["drug_name"]]
              .merge(integrated, on="drug_name", how="left")
              .sort_values("integrated", ascending=False)
              .reset_index(drop=True))
    ranked.to_csv(out_dir / "ranked_drugs.csv", index=False)
    print(f"Ranked candidates: {len(ranked)} drugs -> {out_dir / 'ranked_drugs.csv'}")

    n_wcs = len(pd.read_csv(WORKFLOW_ROOT / "10_filter" / "output" / "wcs_filtered.csv"))
    write_report(ranked, len(integrated), n_wcs, integrated["integrated"].values)

    merge_cols = (["drug_name", "WCS", "RGES", "Zhang_c", "pathway_reversal", "integrated"]
                  + [f"{n}_efficacy" for n in critical])
    merged = ranked[merge_cols]
    render_rank_heatmap(
        merged, integrated, critical,
        f"MIRI Drug Ranking — Score Percentiles ({len(merged)} filtered drugs)",
        "heatmap_rank.png")

    integ_pct = therapeutic_percentile_matrix(merged, integrated, ["integrated"])[:, 0]
    top = merged[integ_pct >= 0.5]
    render_rank_heatmap(
        top, integrated, critical,
        f"MIRI Drug Ranking — Integrated pct >= 0.5 ({len(top)}/{len(merged)} filtered drugs)",
        "heatmap_rank@50.png")

    print(f"\nStep 11 complete — top candidates:")
    for i, (_, row) in enumerate(ranked.head(10).iterrows()):
        print(f"    {i+1:2d}. {row['drug_name']:40s} integrated={row['integrated']:+.4f}")


if __name__ == "__main__":
    main()
