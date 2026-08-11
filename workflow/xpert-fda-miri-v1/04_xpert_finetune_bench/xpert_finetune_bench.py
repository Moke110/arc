#!/usr/bin/env python3
"""
Step 04 — Benchmark XPert on DTOXS splits.

Finetunes XPert on the evaluation DToxS splits (warm-start, cold-drug,
cold-cate; split columns created by create_dtoxs_splits.py in input/dtoxs.h5ad)
and benchmarks each fold on its held-out test set. Idempotent: splits with
existing fold weights are not retrained, and splits with existing per-fold
result CSVs are not re-evaluated (pass --force to redo both).

XPert writes training/test artifacts under <workflow>/finetune/output/ (the
library derives this staging path from the config root); this step moves them
into 04_xpert_finetune_bench/output/ when each run finishes.

Outputs:
  04_xpert_finetune_bench/output/dtoxs_<split>/   — fold weights + per-fold metrics
  04_xpert_finetune_bench/output/benchmark.csv    — per-fold test metrics, all splits
  04_xpert_finetune_bench/benchmark.png           — mean +/- std per split (deg metrics)

Usage:
    python 04_xpert_finetune_bench/xpert_finetune_bench.py            # idempotent
    python 04_xpert_finetune_bench/xpert_finetune_bench.py --force    # retrain + re-evaluate
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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, XPERT_ROOT  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "04_xpert_finetune_bench"
OUTPUT_DIR = STEP_DIR / "output"
CONFIG = WORKFLOW_ROOT / "05_xpert_finetune" / "config" / "dtoxs"
STAGING_ROOT = WORKFLOW_ROOT / "finetune" / "output"
DEFAULT_PRETRAINED = (
    WORKFLOW_ROOT / "03_xpert_pretrain" / "pretrain_mdmt_full_50_epoch.pth"
)
SPLITS = ["warm-start", "cold-drug", "cold-cate"]
PLOT_METRICS = ["MSE_deg", "RMSE_deg", "MAE_deg",
                "Pearson_deg", "Spearman_deg", "R2_deg"]


def run_xpert(mode, split, device, extra=()):
    cmd = [
        sys.executable, "-m", "src.train",
        "--mode", mode,
        "--model", "XPert",
        "--config", str(CONFIG),
        "--drug_feat", "unimol",
        "--drug_embeds", "unimol,HG",
        "--split_method", split,
        "--dataset", "dtoxs",
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
    out_dir = OUTPUT_DIR / f"dtoxs_{split}"
    staging = STAGING_ROOT / f"dtoxs_{split}"
    if not force and list(out_dir.glob("*.pth")):
        print(f"Cached fold weights present: {out_dir} — skipping training")
        return
    shutil.rmtree(staging, ignore_errors=True)
    run_xpert("train", split, device,
              extra=["--pretrained_model", pretrained, "--num_epochs", "4"])
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
    out_dir = OUTPUT_DIR / f"dtoxs_{split}"
    staging = STAGING_ROOT / f"dtoxs_{split}"
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


def aggregate():
    rows = []
    for split in SPLITS:
        out_dir = OUTPUT_DIR / f"dtoxs_{split}"
        for csv in sorted(out_dir.glob("*_train_result.csv")):
            fold = csv.stem.split("_fold-")[-1].replace("_train_result", "")
            df = pd.read_csv(csv)
            row = df.iloc[0].to_dict()  # tr/val/test rows are identical in test mode
            row.update({"split": split, "fold": int(fold)})
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
    x = np.arange(len(SPLITS))
    for ax, metric in zip(axes.flatten(), metrics):
        means = [bench.loc[bench["split"] == s, metric].mean() for s in SPLITS]
        stds = [bench.loc[bench["split"] == s, metric].std() for s in SPLITS]
        ax.bar(x, means, yerr=stds, capsize=4, color="steelblue", alpha=0.85,
               edgecolor="white")
        ax.set_xticks(x)
        ax.set_xticklabels(SPLITS, fontsize=10)
        ax.set_title(metric, fontsize=12, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("XPert DTOXS finetune benchmark — test metrics per split "
                 "(mean ± std across folds)", fontsize=14, fontweight="bold")
    fig.tight_layout()
    out_png = STEP_DIR / "benchmark.png"
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")


def main():
    parser = argparse.ArgumentParser(description="Benchmark XPert on DTOXS splits")
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

    for split in SPLITS:
        train_split(split, args.device, args.pretrained_model, args.force)
        eval_split(split, args.device, args.force or args.force_eval)

    shutil.rmtree(WORKFLOW_ROOT / "finetune", ignore_errors=True)

    bench = aggregate()
    plot_benchmark(bench)
    print(f"Step 04 complete — benchmark in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
