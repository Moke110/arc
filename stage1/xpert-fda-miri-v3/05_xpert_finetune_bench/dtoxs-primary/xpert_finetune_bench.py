#!/usr/bin/env python3
"""
Step 05 — Benchmark XPert on the union CM splits (ADR 0015) — dtoxs-primary.

Finetunes XPert on the evaluation splits — warm-start (drug-overlap) and
cold-drug (drug-disjoint), 10 folds each under the unified fold scheme
(ADR 0014; split columns imported from v1.1 for the primary arm + iPSC
isolation by create_cm_splits.py in the union h5ad). The benchmark trains on
the workflow-built union h5ad (03_cm_finetune_data/output/cm_union.h5ad),
not the DToxS-only frozen h5ad. In addition it evaluates the workflow-owned
L1000 pretrain weight (zero-shot, no finetuning) on every fold test set under
config/cm_pretrain_eval.yaml, and plots that as the FIRST column in
dtoxs-primary-bench.png. Idempotent: folds with existing weights are not
retrained, and folds with existing per-fold result CSVs are not re-evaluated
(pass --force to redo both).

Because every fold's test set comes from the dtoxs_primary arm (the
gse244740/gse262419 rows are train-only, ADR 0017), the benchmark's outputs
live under 05_xpert_finetune_bench/dtoxs-primary/ and its plot is named
dtoxs-primary-bench.png (ADR 0018). The data-ablation experiment
(xpert_data_ablation.py) measures the finetune-corpus scale with its own
decoupled per-seed unified test sets and must never be mixed with these folds.

XPert writes training/test artifacts under <workflow>/finetune/output/ (the
library derives this staging path from the config root); this step moves them
into 05_xpert_finetune_bench/dtoxs-primary/ when each run finishes.

Outputs:
  05_xpert_finetune_bench/dtoxs-primary/cm_warm-start/    — fold weights + per-fold metrics
  05_xpert_finetune_bench/dtoxs-primary/pretrain_warm-start/ — L1000-pretrain zero-shot fold metrics
  05_xpert_finetune_bench/dtoxs-primary/benchmark.csv     — per-fold test metrics, all splits
  05_xpert_finetune_bench/dtoxs-primary/dtoxs-primary-bench.png — pretrain first, then mean +/- std per split

Usage:
    python 05_xpert_finetune_bench/dtoxs-primary/xpert_finetune_bench.py            # idempotent
    python 05_xpert_finetune_bench/dtoxs-primary/xpert_finetune_bench.py --force    # retrain + re-evaluate
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, XPERT_ROOT  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "05_xpert_finetune_bench"
# The official fold benchmark: its fold test sets all come from the
# dtoxs_primary arm, so everything lives under dtoxs-primary/ and the plot is
# dtoxs-primary-bench.png (ADR 0018).
OUTPUT_DIR = STEP_DIR / "dtoxs-primary"
CONFIG = WORKFLOW_ROOT / "06_xpert_finetune" / "config" / "cm"
PRETRAIN_EVAL_CONFIG = STEP_DIR / "config" / "cm_pretrain_eval"
STAGING_ROOT = WORKFLOW_ROOT / "finetune" / "output"
DEFAULT_PRETRAINED = (
    WORKFLOW_ROOT / "04_xpert_pretrain" / "pretrain_mdmt_full_50_epoch.pth"
)
SPLITS = ["warm-start", "cold-drug"]
PRETRAIN_SPLIT = "pretrain"
PLOT_GROUPS = [PRETRAIN_SPLIT] + SPLITS
PLOT_COLORS = {PRETRAIN_SPLIT: "#e8a13c", "warm-start": "steelblue",
               "cold-drug": "steelblue"}
N_FOLDS = 10
PLOT_METRICS = ["MSE_deg", "RMSE_deg", "MAE_deg",
                "Pearson_deg", "Spearman_deg", "R2_deg"]
NUM_EPOCHS = 4


def is_test_csv(csv):
    # Only accept mode=test outputs (uniform test-metric rows). Training-time
    # CSVs carry [train, val, val] rows (src.train overwrites test with val),
    # which would silently corrupt the aggregate — re-evaluate those folds.
    df = pd.read_csv(csv)
    return len(df) >= 3 and df.iloc[0].equals(df.iloc[2])


def run_xpert(mode, split, device, extra=(), config=CONFIG):
    cmd = [
        sys.executable, "-m", "src.train",
        "--mode", mode,
        "--model", "XPert",
        "--config", str(config),
        "--drug_feat", "unimol",
        "--drug_embeds", "unimol,HG",
        "--split_method", split,
        "--dataset", "cm",
        "--use_gradscaler", "True",
        "--include_cell_idx", "True",
        "--batch_size", "32",
        "--init_epoch", "0",
        "--finetune_lr", "0.002",
        "--device", device,
        *extra,
    ]
    print(f"Running {mode} ({split}):\n  {' '.join(cmd)}")
    subprocess.run(cmd, cwd=XPERT_ROOT, check=True)


def train_split(split, device, pretrained, force):
    out_dir = OUTPUT_DIR / f"cm_{split}"
    staging = STAGING_ROOT / f"cm_{split}"
    if not force and list(out_dir.glob("*.pth")):
        print(f"Cached fold weights present: {out_dir} — skipping training")
        return
    shutil.rmtree(staging, ignore_errors=True)
    run_xpert("train", split, device,
              extra=["--pretrained_model", pretrained,
                     "--num_epochs", str(NUM_EPOCHS)])
    if not staging.exists() or not list(staging.glob("*.pth")):
        print(f"ERROR: training produced no weights in {staging}")
        sys.exit(1)
    out_dir.mkdir(parents=True, exist_ok=True)
    for item in staging.iterdir():
        dest = out_dir / item.name
        if dest.exists():
            shutil.rmtree(dest) if dest.is_dir() else dest.unlink()
        shutil.move(str(item), str(dest))
    shutil.rmtree(staging, ignore_errors=True)


def eval_split(split, device, force):
    out_dir = OUTPUT_DIR / f"cm_{split}"
    staging = STAGING_ROOT / f"cm_{split}"
    pths = sorted(out_dir.glob("*_best.pth"))
    if not pths:
        print(f"ERROR: no fold weights in {out_dir} — train first")
        sys.exit(1)
    n_expected = len(pths)
    # Only accept mode=test outputs (uniform test-metric rows). Training-time
    # CSVs carry [train, val, val] rows (src.train overwrites test with val),
    # which would silently corrupt the aggregate — re-evaluate those folds.
    def is_test_csv(csv):
        df = pd.read_csv(csv)
        return len(df) >= 3 and df.iloc[0].equals(df.iloc[2])
    existing = [c for c in out_dir.glob("*_train_result.csv") if is_test_csv(c)]
    if not force and len(existing) >= n_expected:
        print(f"Cached fold metrics present: {out_dir} — skipping evaluation")
        return
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True, exist_ok=True)
    # --pretrained_model is only used by src.train to pick the 'finetune'
    # staging root; the tested weights come from --saved_model_path.
    run_xpert("test", split, device,
              extra=["--saved_model_path", str(out_dir),
                     "--pretrained_model", str(pths[0])])
    for csv in staging.glob("*_train_result.csv"):
        shutil.move(str(csv), str(out_dir / csv.name))
    shutil.rmtree(staging, ignore_errors=True)
    remaining = list(out_dir.glob("*_train_result.csv"))
    if len(remaining) < n_expected:
        print(f"WARNING: expected {n_expected} fold result CSVs, "
              f"got {len(remaining)} in {out_dir}")


def eval_pretrain_split(split, device, pretrained, force):
    """Zero-shot evaluation of the L1000 pretrain weight (no finetuning).

    Runs src.train --mode test under config/cm_pretrain_eval.yaml with the
    SAME pretrain checkpoint used for every fold, so each fold's held-out test
    set is scored by the un-finetuned base model. Result CSVs land in
    output/pretrain_<split>/ and are aggregated into the 'pretrain' group.
    """
    out_dir = OUTPUT_DIR / f"pretrain_{split}"
    staging = STAGING_ROOT / f"cm_{split}"
    out_dir.mkdir(parents=True, exist_ok=True)
    existing = [c for c in out_dir.glob("*_train_result.csv") if is_test_csv(c)]
    if not force and len(existing) >= N_FOLDS:
        print(f"Cached pretrain fold metrics present: {out_dir} — skipping")
        return
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True, exist_ok=True)
    # --pretrained_model keeps the staging root under finetune/output and is
    # not used for weight selection; --saved_model is the pretrain checkpoint.
    run_xpert("test", split, device,
              extra=["--config", str(PRETRAIN_EVAL_CONFIG),
                     "--saved_model", str(pretrained),
                     "--pretrained_model", str(pretrained)])
    for csv in staging.glob("*_train_result.csv"):
        shutil.move(str(csv), str(out_dir / csv.name))
    shutil.rmtree(staging, ignore_errors=True)
    remaining = list(out_dir.glob("*_train_result.csv"))
    if len(remaining) < N_FOLDS:
        print(f"WARNING: expected {N_FOLDS} pretrain fold result CSVs, "
              f"got {len(remaining)} in {out_dir}")


def aggregate():
    rows = []
    for split in SPLITS:
        out_dir = OUTPUT_DIR / f"cm_{split}"
        for csv in sorted(out_dir.glob("*_train_result.csv")):
            fold = csv.stem.split("_fold-")[-1].replace("_train_result", "")
            df = pd.read_csv(csv)
            row = df.iloc[0].to_dict()  # tr/val/test rows are identical in test mode
            row.update({"split": split, "fold": int(fold)})
            rows.append(row)
    # L1000-pretrain zero-shot rows: one 'pretrain' group across all 20 fold
    # test sets (warm-start folds 1-10, cold-drug folds 11-20).
    for split in SPLITS:
        out_dir = OUTPUT_DIR / f"pretrain_{split}"
        for csv in sorted(out_dir.glob("*_train_result.csv")):
            fold = csv.stem.split("_fold-")[-1].replace("_train_result", "")
            df = pd.read_csv(csv)
            row = df.iloc[0].to_dict()
            offset = N_FOLDS if split == "cold-drug" else 0
            row.update({"split": PRETRAIN_SPLIT, "fold": int(fold) + offset})
            rows.append(row)
    if not rows:
        print("ERROR: no per-fold result CSVs found — run evaluation first")
        sys.exit(1)
    bench = pd.DataFrame(rows)
    front = ["split", "fold"]
    bench = bench[front + [c for c in bench.columns if c not in front]]
    bench = bench.sort_values(["split", "fold"]).reset_index(drop=True)
    out_csv = OUTPUT_DIR / "benchmark.csv"
    bench.to_csv(out_csv, index=False)
    print(f"Aggregated {len(bench)} folds -> {out_csv}")
    return bench


def plot_benchmark(bench):
    metrics = [m for m in PLOT_METRICS if m in bench.columns]
    if not metrics:
        print("WARNING: none of the plot metrics found in benchmark.csv — skipping plot")
        return
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    x = np.arange(len(PLOT_GROUPS))
    rng = np.random.default_rng(2024)
    for ax, metric in zip(axes.flatten(), metrics):
        means = [bench.loc[bench["split"] == g, metric].mean() for g in PLOT_GROUPS]
        stds = [bench.loc[bench["split"] == g, metric].std() for g in PLOT_GROUPS]
        ax.bar(x, means, yerr=stds, capsize=4,
               color=[PLOT_COLORS[g] for g in PLOT_GROUPS], alpha=0.85,
               edgecolor="white")
        for g_idx, group in enumerate(PLOT_GROUPS):
            vals = bench.loc[bench["split"] == group, metric].astype(float).values
            ax.scatter(x[g_idx] + rng.uniform(-0.12, 0.12, size=len(vals)),
                       vals, color="black", s=14, zorder=3, alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(PLOT_GROUPS, fontsize=10)
        ax.set_title(metric, fontsize=12, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    legend = [plt.Rectangle((0, 0), 1, 1, fc=PLOT_COLORS[PRETRAIN_SPLIT]),
              plt.Rectangle((0, 0), 1, 1, fc=PLOT_COLORS["warm-start"])]
    fig.legend(legend, [f"{PRETRAIN_SPLIT} (L1000, no finetune)",
                        "union-CM finetune"],
               loc="lower center", ncol=2, fontsize=11,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("XPert union-CM finetune benchmark — test metrics per group "
                 "(mean ± std across folds; pretrain = L1000 weight, "
                 "no finetune)", fontsize=14, fontweight="bold")
    fig.tight_layout()
    out_png = OUTPUT_DIR / "dtoxs-primary-bench.png"
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")


def main():
    parser = argparse.ArgumentParser(description="Benchmark XPert on the union CM splits")
    parser.add_argument("--pretrained-model", default=str(DEFAULT_PRETRAINED),
                        help="Pretrained XPert checkpoint (for training missing folds)")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--force", action="store_true",
                        help="Retrain and re-evaluate even if cached artifacts exist")
    parser.add_argument("--force-eval", action="store_true",
                        help="Re-evaluate existing weights even if cached metrics exist")
    args = parser.parse_args()

    if not Path(str(CONFIG) + ".yaml").exists():
        print(f"ERROR: config not found: {CONFIG}")
        sys.exit(1)
    if not Path(str(PRETRAIN_EVAL_CONFIG) + ".yaml").exists():
        print(f"ERROR: pretrain eval config not found: {PRETRAIN_EVAL_CONFIG}")
        sys.exit(1)

    for split in SPLITS:
        train_split(split, args.device, args.pretrained_model, args.force)
        eval_split(split, args.device, args.force or args.force_eval)

    for split in SPLITS:
        eval_pretrain_split(split, args.device, args.pretrained_model,
                            args.force or args.force_eval)

    shutil.rmtree(WORKFLOW_ROOT / "finetune", ignore_errors=True)

    bench = aggregate()
    plot_benchmark(bench)
    print(f"Step 05 complete — benchmark in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
