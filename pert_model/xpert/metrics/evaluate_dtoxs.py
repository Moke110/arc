#!/usr/bin/env python
"""Evaluate pretrain and finetune XPert models on DTOXS dataset.

Generates per-metric CSV files and PNG bar charts in xpert/metrics/dtoxs/.

Groups:
  L1000              — pretrain_mdmt_full_50_epoch.pth evaluated on 10 warm-start-fold test sets
  DToxS warm-start   — dtoxs_warm-start_fold-*_best.pth (10 folds)
  DToxS cold-drug    — dtoxs_cold-drug_fold-*_best.pth (10 folds)
  DToxS cold-cate    — dtoxs_cold-cate_fold-*_best.pth (6 folds)
"""

import sys
import os
import logging
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import torch
import scanpy as sc
from torch.utils.data import DataLoader
from tqdm import tqdm

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import load_config
from src.dataset import MyDataset
from src.model_xpert import XPertNet
from src.metrics import get_metrics_infer

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# Pretrain weight is workflow-owned since ADR 0008
PRETRAIN_CKPT = PROJECT_ROOT / "../../workflow/xpert-fda-miri-v1/03_xpert_pretrain/pretrain_mdmt_full_50_epoch.pth"
FINETUNE_WARM_START_DIR = PROJECT_ROOT / "finetune/output/dtoxs_warm-start"
FINETUNE_COLD_DRUG_DIR = PROJECT_ROOT / "finetune/output/dtoxs_cold-drug"
FINETUNE_COLD_CATE_DIR = PROJECT_ROOT / "finetune/output/dtoxs_cold-cate"
OUTPUT_DIR = PROJECT_ROOT / "metrics/dtoxs"
H5AD_PATH = PROJECT_ROOT / "finetune/data/dtoxs.h5ad"

# Metric names matching existing PNGs
METRIC_NAMES = [
    "MAE_deg", "MSE_deg", "Pearson_deg",
    "Precision@10 Negative_deg", "Precision@10 Positive_deg",
    "Precision@20 Negative_deg", "Precision@20 Positive_deg",
    "Precision@50 Negative_deg", "Precision@50 Positive_deg",
    "Precision@100 Negative_deg", "Precision@100 Positive_deg",
    "R2_deg", "RMSE_deg", "Spearman_deg",
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def make_args(dataset, mode="infer"):
    """Build args namespace for a given dataset."""
    return argparse.Namespace(
        dataset=dataset,
        drug_feat="unimol",
        drug_embeds="unimol,HG",
        mode=mode,
        pretrained_mode="global",
        include_cell_idx=True,
        output_attention=False,
        output_cls_embed=False,
        wo_HG=False,
        wo_atom=False,
        wo_atom_HG=False,
        wo_unimol=False,
        wo_ppi=False,
        use_gene_pos_emed=False,
        batch_size=256,
    )


def get_test_dataloader(adata, data_args, config, logger, fold_col):
    """Create a DataLoader for the test split of a given fold column."""
    test_data = adata[adata.obs[fold_col] == "test"].copy()
    max_value = config["dataset"]["max_value"]
    min_value = config["dataset"]["min_value"]
    drug_feat = np.load(config["drug_embeddings"]["unimol"], allow_pickle=True)

    dataset = MyDataset(test_data, drug_feat, args=data_args, config=config,
                        logger=logger, max_value=max_value, min_value=min_value)
    return DataLoader(dataset, batch_size=data_args.batch_size, shuffle=False, num_workers=0)


def run_inference(model, dataloader, device):
    """Run inference, return y_true (x_pert), y_pred (x_pert), ctl_true."""
    model.eval()
    y_true_ls, y_pred_ls, ctl_true_ls = [], [], []

    with torch.no_grad():
        for data in tqdm(dataloader, desc="Infer", leave=False):
            trt_output, ctl_output, deg_output, trt_raw, ctl_raw, _, _, _, _ = model(data)
            y_true_ls.append(trt_raw.cpu().numpy())
            y_pred_ls.append((deg_output.cpu() + ctl_raw.cpu()).numpy())
            ctl_true_ls.append(ctl_raw.cpu().numpy())

    y_true = np.concatenate(y_true_ls, axis=0)
    y_pred = np.concatenate(y_pred_ls, axis=0)
    ctl_true = np.concatenate(ctl_true_ls, axis=0)
    return y_true, y_pred, ctl_true


def build_model(args, config, device, logger, ckpt_path, strict=True):
    """Load model from checkpoint."""
    ckpt = torch.load(ckpt_path, map_location=device, weights_only=False)
    model = XPertNet(args, config, device, logger).to(device)
    model.load_state_dict(ckpt["model_state_dict"], strict=strict)
    model.args.mode = "infer"
    return model


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def evaluate_models(output_dir="metrics/dtoxs"):
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    OUTPUT_DIR = PROJECT_ROOT / output_dir
    print(f"Using device: {device}")

    # Two configs: pretrain (for pretrain model architecture) and dtoxs (for data + finetuned models)
    dtoxs_config = load_config("finetune/config/dtoxs")
    pretrain_config = load_config("pretrain/config/l1000")

    # Two arg sets: pretrain model needs 'mdmt' in dataset name; finetune needs 'dtoxs'
    pretrain_args = make_args("l1000_mdmt_full")
    dtoxs_args = make_args("dtoxs")
    # DataLoader needs mode='train' (not 'infer') for 10-item tuple instead of 14
    loader_args = make_args("dtoxs", mode="train")

    # Logger
    logger = logging.getLogger("evaluate")
    logger.setLevel(logging.WARNING)
    if not logger.handlers:
        logger.addHandler(logging.StreamHandler())

    # Load data
    print(f"Loading {H5AD_PATH}...")
    adata = sc.read_h5ad(H5AD_PATH)

    # Container: {metric_name: {group_name: [fold_values]}}
    all_results: dict = {
        m: {"L1000": [], "DToxS warm-start": [], "DToxS cold-drug": [], "DToxS cold-cate": []}
        for m in METRIC_NAMES
    }

    # -----------------------------------------------------------------------
    # 1. Pretrain model on 10 warm-start folds
    # -----------------------------------------------------------------------
    print("\n=== PRETRAIN model on warm-start folds ===")
    model = build_model(pretrain_args, pretrain_config, device, logger, str(PRETRAIN_CKPT), strict=False)

    for fold_i in range(1, 11):
        fold_col = f"warm-start_fold-{fold_i}"
        print(f"  Fold {fold_i}/10 ({fold_col})  test N={len(adata[adata.obs[fold_col] == 'test'])}")
        loader = get_test_dataloader(adata, loader_args, dtoxs_config, logger, fold_col)
        y_true, y_pred, ctl_true = run_inference(model, loader, device)
        metrics, _ = get_metrics_infer(y_true, y_pred, ctl_true)
        for m in METRIC_NAMES:
            all_results[m]["L1000"].append(metrics.get(m, np.nan))

    del model

    # -----------------------------------------------------------------------
    # 2. Finetuned warm-start models (10 folds)
    # -----------------------------------------------------------------------
    print("\n=== FINETUNE warm-start models ===")
    for fold_i in range(1, 11):
        ckpt_path = FINETUNE_WARM_START_DIR / f"dtoxs_warm-start_fold-{fold_i}_best.pth"
        print(f"  Fold {fold_i}/10 ({ckpt_path.name})  test N={len(adata[adata.obs[f'warm-start_fold-{fold_i}'] == 'test'])}")
        model = build_model(dtoxs_args, dtoxs_config, device, logger, str(ckpt_path))

        loader = get_test_dataloader(adata, loader_args, dtoxs_config, logger, f"warm-start_fold-{fold_i}")
        y_true, y_pred, ctl_true = run_inference(model, loader, device)
        metrics, _ = get_metrics_infer(y_true, y_pred, ctl_true)
        for m in METRIC_NAMES:
            all_results[m]["DToxS warm-start"].append(metrics.get(m, np.nan))
        del model

    # -----------------------------------------------------------------------
    # 3. Finetuned cold-drug models (10 folds)
    # -----------------------------------------------------------------------
    print("\n=== FINETUNE cold-drug models ===")
    for fold_i in range(1, 11):
        ckpt_path = FINETUNE_COLD_DRUG_DIR / f"dtoxs_cold-drug_fold-{fold_i}_best.pth"
        fold_col = f"cold-drug_fold-{fold_i}"
        print(f"  Fold {fold_i}/10 ({ckpt_path.name})  test N={len(adata[adata.obs[fold_col] == 'test'])}")
        model = build_model(dtoxs_args, dtoxs_config, device, logger, str(ckpt_path))

        loader = get_test_dataloader(adata, loader_args, dtoxs_config, logger, fold_col)
        y_true, y_pred, ctl_true = run_inference(model, loader, device)
        metrics, _ = get_metrics_infer(y_true, y_pred, ctl_true)
        for m in METRIC_NAMES:
            all_results[m]["DToxS cold-drug"].append(metrics.get(m, np.nan))
        del model

    # -----------------------------------------------------------------------
    # 4. Finetuned cold-cate models (6 folds)
    # -----------------------------------------------------------------------
    print("\n=== FINETUNE cold-cate models ===")
    for fold_i in range(1, 7):
        ckpt_path = FINETUNE_COLD_CATE_DIR / f"dtoxs_cold-cate_fold-{fold_i}_best.pth"
        fold_col = f"cold-cate_fold-{fold_i}"
        print(f"  Fold {fold_i}/6 ({ckpt_path.name})  test N={len(adata[adata.obs[fold_col] == 'test'])}")
        model = build_model(dtoxs_args, dtoxs_config, device, logger, str(ckpt_path))

        loader = get_test_dataloader(adata, loader_args, dtoxs_config, logger, fold_col)
        y_true, y_pred, ctl_true = run_inference(model, loader, device)
        metrics, _ = get_metrics_infer(y_true, y_pred, ctl_true)
        for m in METRIC_NAMES:
            all_results[m]["DToxS cold-cate"].append(metrics.get(m, np.nan))
        del model

    # -----------------------------------------------------------------------
    # Save CSVs
    # -----------------------------------------------------------------------
    print("\n=== Saving CSVs ===")
    group_order = ["L1000", "DToxS warm-start", "DToxS cold-drug", "DToxS cold-cate"]
    max_rows = max(len(all_results[m][g]) for m in METRIC_NAMES for g in group_order)

    for metric in METRIC_NAMES:
        csv_path = OUTPUT_DIR / f"{metric}.csv"
        df_data = {}
        for group in group_order:
            vals = all_results[metric][group]
            padded = list(vals) + [np.nan] * (max_rows - len(vals))
            df_data[group] = padded
        df = pd.DataFrame(df_data)
        df.to_csv(csv_path, index=False)
        print(f"  {csv_path}")

    # -----------------------------------------------------------------------
    # Generate PNG bar charts
    # -----------------------------------------------------------------------
    print("\n=== Generating PNGs ===")
    plt.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "figure.dpi": 100,
    })

    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    for metric in METRIC_NAMES:
        fig, ax = plt.subplots(figsize=(8, 6))

        group_means = []
        group_stds = []
        group_labels = []

        for gi, group in enumerate(group_order):
            vals = np.array(all_results[metric][group])
            vals = vals[~np.isnan(vals)]
            group_means.append(np.mean(vals))
            group_stds.append(np.std(vals, ddof=1) if len(vals) > 1 else 0)
            group_labels.append(group)

        x = np.arange(len(group_labels))
        bars = ax.bar(x, group_means, yerr=group_stds, color=colors, capsize=6,
                       edgecolor="white", linewidth=0.8)

        for bar, mean_val in zip(bars, group_means):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.003,
                    f"{mean_val:.3f}", ha="center", va="bottom", fontsize=9)

        ax.set_xticks(x)
        ax.set_xticklabels(group_labels, rotation=15, ha="right")
        ax.set_ylabel(metric.replace("_", " "))
        ax.set_title(metric.replace("_", " "))
        ax.grid(axis="y", alpha=0.3)

        fig.tight_layout()
        png_path = OUTPUT_DIR / f"{metric}.png"
        fig.savefig(png_path, dpi=100, bbox_inches="tight")
        plt.close(fig)
        print(f"  {png_path}")

    print("\nDone.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate XPert models on DTOXS dataset")
    parser.add_argument("--output-dir", type=str, default="metrics/dtoxs",
                        help="Output directory for metrics (default: metrics/dtoxs)")
    args = parser.parse_args()
    evaluate_models(args.output_dir)
