#!/usr/bin/env python
"""Generate warm-start and cm_full split columns in the union h5ad.

Creates two new sets of .obs columns on 03_cm_finetune_data/output/cm_union.h5ad:

  warm-start (10 folds, evaluation-only):
    - Drug-stratified: each drug's samples split ~80/20 train/test per fold
    - No validation set
    - Each sample appears in the test set of >= 1 fold
    - Columns: warm-start_fold-1 through warm-start_fold-10
    - Values: 'train' or 'test'

  cm_full (1 fold, prediction-only):
    - All samples assigned 'train' (no test, no valid)
    - Column: cm_full
    - Values: 'train'

Only samples with a valid HG-resolved pert_idx (0 <= pert_idx < N_L1000_DRUGS) are
included in splits; non-small-molecule / excluded perturbations are left
empty ('').

Usage:
    python 05_xpert_finetune_bench/create_cm_splits.py
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "_lib"))
from utils import N_L1000_DRUGS  # noqa: E402

logger = logging.getLogger("create_cm_splits")
logger.setLevel(logging.INFO)
if not logger.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)


def create_warm_start_splits(adata, n_folds=10, test_frac=0.2, seed=2024):
    """Create warm-start split columns: drug-stratified, no validation, each sample in >=1 test fold.

    For each fold, ~test_frac of each drug's samples are randomly assigned to 'test',
    the remaining ~(1 - test_frac) to 'train'. No 'valid' split.

    Args:
        adata: AnnData with .obs containing 'pert_idx' column
        n_folds: number of folds (default 10)
        test_frac: fraction of each drug's samples to assign to test (default 0.2)
        seed: random seed for reproducibility

    Returns:
        adata with new .obs columns warm-start_fold-1 through warm-start_fold-N
    """
    rng = np.random.default_rng(seed)
    pert_idx = adata.obs['pert_idx'].values.astype(int)

    # Only assign splits for valid drug indices (small molecules in HG range)
    valid_mask = (pert_idx >= 0) & (pert_idx < N_L1000_DRUGS)
    valid_indices = np.where(valid_mask)[0]

    drug_to_samples = {}
    for i in valid_indices:
        did = pert_idx[i]
        drug_to_samples.setdefault(did, []).append(i)

    drug_ids = sorted(drug_to_samples.keys())
    logger.info(f"Found {len(drug_ids)} unique drugs with {len(valid_indices)} total samples")
    logger.info(f"Total samples in adata: {adata.n_obs} "
                f"({np.sum(~valid_mask)} non-HG excluded)")

    sample_test_count = np.zeros(adata.n_obs, dtype=int)

    for fold_i in range(1, n_folds + 1):
        col_name = f"warm-start_fold-{fold_i}"
        # dtype must fit 'train'/'test' — np.full(..., dtype=str) would give
        # <U1 and truncate to 't'
        fold_values = np.full(adata.n_obs, "", dtype="<U10")

        for did, sample_indices in drug_to_samples.items():
            n_samples = len(sample_indices)
            n_test = max(1, int(np.round(n_samples * test_frac)))
            n_test = min(n_test, n_samples - 1) if n_samples > 1 else 1

            shuffled = rng.permutation(sample_indices)
            test_idx = shuffled[:n_test]
            train_idx = shuffled[n_test:]

            for idx in test_idx:
                fold_values[idx] = "test"
                sample_test_count[idx] += 1
            for idx in train_idx:
                fold_values[idx] = "train"

        adata.obs[col_name] = fold_values
        n_train = np.sum(fold_values == "train")
        n_test = np.sum(fold_values == "test")
        logger.info(f"  {col_name}: train={n_train}, test={n_test}")

    uncovered = []
    for i in valid_indices:
        if sample_test_count[i] == 0:
            uncovered.append(i)
    if uncovered:
        logger.warning(f"{len(uncovered)} samples never assigned to test! Fixing...")
        for i in uncovered:
            rescue_fold = rng.integers(1, n_folds + 1)
            adata.obs.loc[i, f"warm-start_fold-{rescue_fold}"] = "test"
            sample_test_count[i] += 1

    coverage = sample_test_count[valid_indices]
    logger.info(f"Test fold coverage: min={coverage.min()}, max={coverage.max()}, "
                f"mean={coverage.mean():.1f}, median={np.median(coverage):.0f}")
    return adata


def create_cm_full_split(adata):
    """Create cm_full column: all valid samples = 'train', others = ''."""
    pert_idx = adata.obs['pert_idx'].values.astype(int)
    fold_values = np.where((pert_idx >= 0) & (pert_idx < N_L1000_DRUGS), "train", "")
    adata.obs["cm_full"] = fold_values
    n_train = np.sum(fold_values == "train")
    logger.info(f"  cm_full: train={n_train}, test=0, valid=0")
    return adata


def main():
    parser = argparse.ArgumentParser(description="Generate warm-start and cm_full splits in cm_union.h5ad")
    parser.add_argument("--data", default="03_cm_finetune_data/output/cm_union.h5ad",
                        help="Path to the union h5ad relative to the workflow root")
    parser.add_argument("--seed", type=int, default=2024, help="Random seed")
    parser.add_argument("--n_folds", type=int, default=10, help="Number of warm-start folds")
    args = parser.parse_args()

    data_path = PROJECT_ROOT / args.data
    if not data_path.exists():
        logger.error(f"Data file not found: {data_path}")
        sys.exit(1)

    logger.info(f"Loading {data_path}...")
    adata = sc.read_h5ad(data_path)
    logger.info(f"Loaded: {adata.n_obs} obs, {adata.n_vars} vars")

    # Normalize the obs index to a clean int RangeIndex so the split rescue
    # loop can use .loc[i] without growing the frame (the union h5ad is written
    # with a string index).
    import pandas as pd
    adata.obs.index = pd.RangeIndex(adata.n_obs)

    if 'pert_idx' not in adata.obs.columns:
        logger.error("'pert_idx' column not found in adata.obs — cannot create splits")
        sys.exit(1)

    logger.info("\n=== Creating warm-start splits (10 folds, drug-stratified, no validation) ===")
    adata = create_warm_start_splits(adata, n_folds=args.n_folds, seed=args.seed)

    logger.info("\n=== Creating cm_full split (all-train, prediction only) ===")
    adata = create_cm_full_split(adata)

    logger.info(f"\nSaving to {data_path}...")
    adata.write_h5ad(data_path)
    logger.info("Done. New columns written:")
    for i in range(1, args.n_folds + 1):
        logger.info(f"  warm-start_fold-{i}")
    logger.info("  cm_full")


if __name__ == "__main__":
    main()
