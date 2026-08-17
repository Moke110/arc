#!/usr/bin/env python3
"""
Apples-to-apples benchmark: union-CM warm-start models on v1's DToxS folds.

The union warm-start fold weights (05_xpert_finetune_bench/output/cm_warm-start/
cm_warm-start_fold-{i}_best.pth) are evaluated with mode=test on v1's frozen
finetune h5ad (workflow/xpert-fda-miri-v1/input/dtoxs.h5ad) using its own
warm-start_fold-1..10 test masks. After the union data-alignment (author time
map {3,6,24}, 24 h -> time_idx 2, quantile alignment to the pretrain range
[0, 15.002]), the union's semantics match v1's DToxS encoding, so the two
models are compared on the IDENTICAL test set — the MSE_deg / Pearson_deg /
... gap is then purely the union finetune's effect.

Test-only config: 05_xpert_finetune_bench/config/cm_eval_dtoxs.yaml
  (dataset name `cm`, DTOX_data_root -> v1 dtoxs.h5ad, cm block, [0, 15.002]
  bins, warm-start folds = v1's warm-start_fold-1..10 columns).

Outputs:
  05_xpert_finetune_bench/output/cm_eval_dtoxs_warm-start/ — per-fold test CSVs
  05_xpert_finetune_bench/output/benchmark_dtoxs_folds_vs_v1.csv
      per-fold table + mean ± std vs v1's archived warm-start benchmark

Usage:
    python 05_xpert_finetune_bench/eval_dtoxs_folds.py            # idempotent
    python 05_xpert_finetune_bench/eval_dtoxs_folds.py --force    # re-evaluate
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, V1_WORKFLOW_ROOT, XPERT_ROOT  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "05_xpert_finetune_bench"
OUTPUT_DIR = STEP_DIR / "output"
EVAL_OUT_DIR = OUTPUT_DIR / "cm_eval_dtoxs_warm-start"
CONFIG = STEP_DIR / "config" / "cm_eval_dtoxs"
WEIGHTS_DIR = OUTPUT_DIR / "cm_warm-start"
STAGING_ROOT = WORKFLOW_ROOT / "finetune" / "output" / "cm_warm-start"
METRICS = ["MSE_deg", "RMSE_deg", "MAE_deg",
           "Pearson_deg", "Spearman_deg", "R2_deg"]


def run_eval(device, force):
    pths = sorted(WEIGHTS_DIR.glob("cm_warm-start_fold-*_best.pth"))
    if not pths:
        print(f"ERROR: no union warm-start fold weights in {WEIGHTS_DIR} — "
              f"run 05 xpert_finetune_bench first")
        sys.exit(1)
    if not force and list(EVAL_OUT_DIR.glob("*_train_result.csv")):
        print(f"Cached eval metrics present: {EVAL_OUT_DIR} — skipping "
              f"(pass --force to re-evaluate)")
        return
    shutil.rmtree(STAGING_ROOT, ignore_errors=True)
    STAGING_ROOT.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-m", "src.train",
        "--mode", "test",
        "--model", "XPert",
        "--config", str(CONFIG),
        "--drug_feat", "unimol",
        "--drug_embeds", "unimol,HG",
        "--split_method", "warm-start",
        "--dataset", "cm",
        "--use_gradscaler", "True",
        "--include_cell_idx", "True",
        "--batch_size", "32",
        "--device", device,
        "--saved_model_path", str(WEIGHTS_DIR),
        "--pretrained_model", str(pths[0]),
    ]
    print(f"Running apples-to-apples eval (union warm-start on v1 DToxS folds):\n"
          f"  {' '.join(cmd)}")
    subprocess.run(cmd, cwd=XPERT_ROOT, check=True)
    EVAL_OUT_DIR.mkdir(parents=True, exist_ok=True)
    for csv in STAGING_ROOT.glob("*_train_result.csv"):
        shutil.move(str(csv), str(EVAL_OUT_DIR / csv.name))
    shutil.rmtree(STAGING_ROOT, ignore_errors=True)
    if not list(EVAL_OUT_DIR.glob("*_train_result.csv")):
        print(f"ERROR: no eval result CSVs produced in {EVAL_OUT_DIR}")
        sys.exit(1)


def aggregate():
    rows = []
    for csv in sorted(EVAL_OUT_DIR.glob("*_train_result.csv")):
        fold = csv.stem.split("_fold-")[-1].replace("_train_result", "")
        df = pd.read_csv(csv)
        row = df.iloc[0].to_dict()  # tr/val/test rows identical in test mode
        row["fold"] = int(fold)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("fold").reset_index(drop=True)


def compare_v1():
    v1 = pd.read_csv(V1_WORKFLOW_ROOT / "04_xpert_finetune_bench"
                     / "output" / "benchmark.csv")
    v1 = v1[v1["split"] == "warm-start"][["fold"] + METRICS] \
        .rename(columns={m: f"v1_{m}" for m in METRICS})
    union = aggregate()
    union = union[["fold"] + METRICS].rename(columns={m: f"v2_{m}" for m in METRICS})
    out = v1.merge(union, on="fold").sort_values("fold")
    out_csv = OUTPUT_DIR / "benchmark_dtoxs_folds_vs_v1.csv"
    out.to_csv(out_csv, index=False)
    print(f"\nPer-fold comparison -> {out_csv}")
    print(out.to_string(index=False))

    rows = []
    for m in METRICS:
        rows.append({
            "metric": m,
            "v2_union_on_dtoxs_mean": out[f"v2_{m}"].mean(),
            "v2_union_on_dtoxs_std": out[f"v2_{m}"].std(),
            "v1_dtoxs_mean": out[f"v1_{m}"].mean(),
            "v1_dtoxs_std": out[f"v1_{m}"].std(),
            "delta_v2_minus_v1": out[f"v2_{m}"].mean() - out[f"v1_{m}"].mean(),
        })
    summary = pd.DataFrame(rows)
    print("\nmean ± std (identical DToxS warm-start test folds):")
    print(summary.to_string(index=False))
    summary_csv = OUTPUT_DIR / "benchmark_dtoxs_folds_summary.csv"
    summary.to_csv(summary_csv, index=False)
    print(f"\nSummary -> {summary_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate union warm-start fold models on v1's DToxS folds")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--force", action="store_true",
                        help="Re-evaluate even if cached eval metrics exist")
    args = parser.parse_args()
    run_eval(args.device, args.force)
    compare_v1()


if __name__ == "__main__":
    main()
