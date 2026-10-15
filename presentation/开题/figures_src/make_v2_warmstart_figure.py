#!/usr/bin/env python3
"""Draw the 开题 slide-24 chart: our own v2 (DToxS-CM) finetuned-model benchmark.

The 【补】 note on 开题3r slide 24 asks for the warm-start PCC, Spearman,
Pos P@20 and Neg P@20 of *our* benchmark — the union-CM finetune benchmark run
by `workflow/xpert-fda-miri-v2/05_xpert_finetune_bench/full/`.

Data: workflow/xpert-fda-miri-v2/05_xpert_finetune_bench/full/output/benchmark.csv
  split == "warm-start"                       (10 drug-stratified folds)
  Pearson_deg                   -> PCC
  Spearman_deg                  -> Spearman
  Precision@20 Positive_deg     -> Pos P@20
  Precision@20 Negative_deg     -> Neg P@20

These `_deg` metrics compare the model's *directly predicted* Delta expression
(`deg_output`; y_pred - ctl) against the true Delta, on the finetune data's own
scale (log2(CPM+1), per-(dataset, plate) quantile-aligned to the L1000 pretrain
reference) — no z-scoring anywhere (see pert_model/xpert/src/train.py
validate(flag=True) and src/metrics.py get_metrics_new). Z-scoring happens only
downstream, on the drug x gene DEG matrix, in 07_xpert_predict/zscore.py.

Bars = mean over the 10 folds, error bars = s.d., black dots = per-fold values
(same convention as the deck's other benchmark figures). Single colour; 1:1.

Output (PNG, 200 dpi) -> presentation/figures/v2_dtoxS_warmstart.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/home/chang/arc")
BENCH = (ROOT / "workflow/xpert-fda-miri-v2/05_xpert_finetune_bench"
         / "full/output/benchmark.csv")
OUT = ROOT / "presentation/figures"
OUT.mkdir(parents=True, exist_ok=True)

SPLIT = "warm-start"
BAR = "#2E6DA4"
METRICS = [
    ("PCC\n($\\Delta$)", "Pearson_deg"),
    ("Spearman", "Spearman_deg"),
    ("Pos P@20", "Precision@20 Positive_deg"),
    ("Neg P@20", "Precision@20 Negative_deg"),
]

plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200,
    "font.size": 13, "axes.labelsize": 15, "axes.titlesize": 13,
    "xtick.labelsize": 13, "ytick.labelsize": 12,
    "axes.linewidth": 1.5, "font.family": "DejaVu Sans",
})


def main():
    df = pd.read_csv(BENCH)
    sub = df[df["split"] == SPLIT]

    labels, means, stds, per_fold = [], [], [], []
    for lab, col in METRICS:
        v = sub[col].astype(float).values
        labels.append(lab); means.append(v.mean()); stds.append(v.std(ddof=1))
        per_fold.append(v)
        print(f"{col:32s} n={len(v)}  mean={v.mean():.3f}  sd={v.std(ddof=1):.3f}")

    fig, ax = plt.subplots(figsize=(4.6, 4.6))
    x = np.arange(len(labels))
    ax.bar(x, means, yerr=stds, capsize=5, width=0.60, color=BAR,
           alpha=0.9, edgecolor="white", linewidth=1.4, zorder=2,
           error_kw=dict(ecolor="#333333", lw=1.5, zorder=4))

    rng = np.random.default_rng(2024)
    for xi, vals in zip(x, per_fold):
        ax.scatter(xi + rng.uniform(-0.13, 0.13, size=len(vals)), vals,
                   color="black", s=22, alpha=0.75, zorder=3)

    for xi, m, s in zip(x, means, stds):
        ax.text(xi, m + s + 0.018, f"{m:.3f}\n±{s:.3f}", ha="center",
                va="bottom", fontsize=11.5, linespacing=0.95)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11.5)
    ax.set_ylabel("score")
    ax.set_ylim(0, 0.60)
    ax.set_yticks(np.arange(0, 0.61, 0.1))
    ax.set_title(f"$\\Delta$ expression · {SPLIT} (n = 10 folds)", pad=8,
                 fontsize=12.5)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    out = OUT / "v2_dtoxS_warmstart.png"
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    print("Saved:", out.relative_to(ROOT))


if __name__ == "__main__":
    main()
