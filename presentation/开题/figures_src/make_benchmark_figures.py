#!/usr/bin/env python3
"""Redraw the slide-14 benchmark figures (开题 3.1) at ~2.5x font.

Two source figures (both from the v2 workflow, which is what the 开题 deck
uses):
  A) 05_xpert_finetune_bench/full/output/benchmark.csv
     "XPert union-CM finetune benchmark — test metrics per group"
     -> bar (mean +/- std) + per-fold scatter, groups pretrain/warm-start/cold-drug
  B) 05_xpert_finetune_bench/data-ablation/ablation.csv
     "finetune-corpus data-ablation — test metrics vs train-size fraction"
     -> per-seed scatter + mean +/- std line vs fraction, 2 modes x 2 metric blocks

Each subplot is redrawn on its own (plus a re-rendered full composite), with
the in-figure text scaled to ~2.5x and the subplot titles unified to the metric
column name (e.g. "Pearson_deg", "JAC").

Outputs (PNG, 200 dpi) -> presentation/figures/
    bench_groups_2.5x.png, bench_group_<metric>.png
    ablation_2.5x.png,     ablation_<mode>_<metric>.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/home/chang/arc")
V2 = ROOT / "workflow/xpert-fda-miri-v2/05_xpert_finetune_bench"
BENCH = V2 / "full/output/benchmark.csv"
ABLATION = V2 / "data-ablation/ablation.csv"
OUT = ROOT / "presentation/figures"
OUT.mkdir(parents=True, exist_ok=True)

BENCH_METRICS = ["MSE_deg", "RMSE_deg", "MAE_deg",
                 "Pearson_deg", "Spearman_deg", "R2_deg"]
GROUPS = ["pretrain", "warm-start", "cold-drug"]
COLORS = {"pretrain": "#e8a13c", "warm-start": "steelblue", "cold-drug": "steelblue"}

ABL_DEG = ["RMSE_deg", "Pearson_deg", "Spearman_deg", "R2_deg"]
ABL_VCC = ["PDS", "JAC", "FID", "REACH"]
ABL_MODES = ["warm-start", "cold-drug"]
FRACTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
X_LABELS = ["0", "20%", "40%", "60%", "80%", "100%"]
ANCHOR_X = min(FRACTIONS) / 2

# ~2.5x the original in-figure text (originals used ~8-14 pt).
plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200,
    "font.size": 28, "axes.labelsize": 30, "axes.titlesize": 30,
    "xtick.labelsize": 26, "ytick.labelsize": 26, "legend.fontsize": 24,
    "axes.titleweight": "bold", "axes.linewidth": 1.8,
    "font.family": "DejaVu Sans",
})


def save(fig, name):
    out = OUT / name
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Saved:", out.relative_to(ROOT))


def bench_panel(ax, bench, metric):
    x = np.arange(len(GROUPS))
    rng = np.random.default_rng(2024)
    means = [bench.loc[bench["split"] == g, metric].mean() for g in GROUPS]
    stds = [bench.loc[bench["split"] == g, metric].std() for g in GROUPS]
    ax.bar(x, means, yerr=stds, capsize=6,
           color=[COLORS[g] for g in GROUPS], alpha=0.85, edgecolor="white")
    for gi, g in enumerate(GROUPS):
        vals = bench.loc[bench["split"] == g, metric].astype(float).values
        ax.scatter(x[gi] + rng.uniform(-0.12, 0.12, size=len(vals)),
                   vals, color="black", s=45, zorder=3, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(GROUPS)
    ax.set_title(metric)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def ablation_panel(ax, abl, mode, metric):
    sub = abl[abl["mode"] == mode]
    anchors = abl[abl["mode"] == "pretrain"]
    x_plot = np.arange(len(X_LABELS))
    means, stds = [], []
    for j, xp in enumerate(x_plot):
        vals = (anchors[metric] if j == 0
                else sub.loc[sub["fraction"] == FRACTIONS[j - 1], metric]
                ).astype(float).values
        means.append(vals.mean()); stds.append(vals.std())
        ax.scatter([xp] * len(vals), vals, s=45, color="black", zorder=3,
                   alpha=0.75)
    means, stds = np.array(means), np.array(stds)
    ax.plot(x_plot, means, "-o", color="steelblue", lw=2.2, zorder=2)
    ax.fill_between(x_plot, means - stds, means + stds, alpha=0.22,
                    color="steelblue")
    ax.set_xticks(x_plot)
    ax.set_xticklabels(X_LABELS)
    ax.tick_params(axis="x", labelsize=20)
    ax.set_title(metric)
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    bench = pd.read_csv(BENCH)
    abl = pd.read_csv(ABLATION)

    # ---- A. benchmark: full composite ----
    fig, axes = plt.subplots(2, 3, figsize=(24, 14))
    for ax, metric in zip(axes.flatten(), BENCH_METRICS):
        bench_panel(ax, bench, metric)
    fig.tight_layout()
    save(fig, "bench_groups_2.5x.png")
    # individual panels
    for metric in BENCH_METRICS:
        fig, ax = plt.subplots(figsize=(9, 8))
        bench_panel(ax, bench, metric)
        save(fig, f"bench_group_{metric}.png")

    # ---- B. data-ablation: full composite ----
    rows = [(m, blk) for m in ABL_MODES for blk in (ABL_DEG, ABL_VCC)]
    fig, axes = plt.subplots(4, 4, figsize=(28, 20))
    for r, (mode, metrics) in enumerate(rows):
        for c, metric in enumerate(metrics):
            ablation_panel(axes[r, c], abl, mode, metric)
        axes[r, 0].set_ylabel(mode, fontweight="bold")
    fig.tight_layout()
    save(fig, "ablation_2.5x.png")
    for mode in ABL_MODES:
        for metric in ABL_DEG + ABL_VCC:
            fig, ax = plt.subplots(figsize=(9, 8))
            ablation_panel(ax, abl, mode, metric)
            ax.set_ylabel(mode)
            save(fig, f"ablation_{mode}_{metric}.png")


if __name__ == "__main__":
    main()
