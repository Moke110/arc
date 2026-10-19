#!/usr/bin/env python3
"""
Step 09 — Filter.

Applies the screening filters to the efficacy components:

  1. WCS filter:      keep drugs with WCS < 0
  2. Pathway filter:  keep drugs with >= MIN_CRITICAL critical pathways whose
                      efficacy > 0

The cardiotoxicity stage is deferred (see docs/adr/0004); no cardiotoxicity
filter is applied. Ranking and integration happen in step 10.

Input:  09_miri_efficacy_scores/output/scores.csv, input/MIRI_pathways.json
Output:
  10_filter/output/wcs_filtered.csv     — WCS < 0 drugs
  10_filter/output/path_filtered.csv    — pathway-filtered drugs
  10_filter/pyramid-filters.png         — filter progression funnel (3 stages)

Usage:
    python 10_filter/filter.py
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT  # noqa: E402
from scoring import load_critical_names  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "10_filter"
SCORES = WORKFLOW_ROOT / "09_miri_efficacy_scores" / "output" / "scores.csv"
PATHWAYS = WORKFLOW_ROOT / "input" / "MIRI_pathways.json"
MIN_CRITICAL = 6


def plot_pyramid(n_all, n_wcs, n_path):
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
                va="center", fontsize=14, fontweight="bold", color="#333333")
    ax.set_yticks(range(len(stages)))
    ax.set_yticklabels([s[0] for s in stages], fontsize=13)
    ax.set_xlim(0, max_w * 1.2)
    ax.invert_yaxis()
    ax.set_xlabel("Number of Drugs", fontsize=12)
    ax.set_title("MIRI Drug Screening: Filter Progression", fontsize=16, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = STEP_DIR / "pyramid-filters.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.relative_to(WORKFLOW_ROOT)}")


def main():
    parser = argparse.ArgumentParser(description="Filter MIRI drug candidates")
    parser.add_argument("--scores", default=str(SCORES))
    parser.add_argument("--pathways", default=str(PATHWAYS))
    parser.add_argument("--min-critical", type=int, default=MIN_CRITICAL,
                        help="Minimum critical pathways with efficacy > 0")
    args = parser.parse_args()

    sns.set_theme(style="whitegrid", context="notebook", palette="muted")

    scores = pd.read_csv(args.scores)
    critical = load_critical_names(args.pathways)
    print(f"Critical pathways ({len(critical)}):")

    # Step 1: WCS < 0
    wcs_pass = scores[scores["WCS"] < 0].copy()
    print(f"  WCS < 0: {len(wcs_pass)}/{len(scores)} drugs")

    # Step 2: pathway filter
    n_positive = pd.Series(0, index=wcs_pass.index)
    for name in critical:
        col = f"{name}_efficacy"
        if col not in wcs_pass.columns:
            print(f"  WARNING: missing column {col} — skipping pathway filter")
            continue
        n_positive += (wcs_pass[col] > 0).astype(int)
    path_pass = wcs_pass[n_positive >= args.min_critical].copy()
    print(f"  >= {args.min_critical}/{len(critical)} critical pathways > 0: "
          f"{len(path_pass)}/{len(wcs_pass)} drugs")

    out_dir = STEP_DIR / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    wcs_pass[["drug_name", "WCS"]].sort_values("WCS").to_csv(
        out_dir / "wcs_filtered.csv", index=False)
    path_pass.assign(n_critical_positive=n_positive[path_pass.index].values)[
        ["drug_name", "n_critical_positive"]].to_csv(
        out_dir / "path_filtered.csv", index=False)

    plot_pyramid(len(scores), len(wcs_pass), len(path_pass))

    print(f"\nStep 09 complete — {len(path_pass)} drugs pass both filters")


if __name__ == "__main__":
    main()
