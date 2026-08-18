#!/usr/bin/env python
"""Generate warm-start and dtoxs split columns in dtoxs.h5ad.

Creates two new sets of .obs columns:

  warm-start (10 folds, evaluation-only):
    - Drug-stratified: each drug's samples split ~80/20 train/test per fold
    - No validation set
    - Each sample appears in the test set of >= 1 fold
    - Columns: warm-start_fold-1 through warm-start_fold-10
    - Values: 'train' or 'test'

  dtoxs (1 fold, prediction-only):
    - All samples assigned 'train' (no test, no valid)
    - Column: dtoxs_full
    - Values: 'train'

Only small-molecule drugs with valid pert_idx (>= 0) are included in splits.
Non-small-molecule perturbations (antibodies, proteins, unmatched) are excluded.

Usage:
    python 05_xpert_finetune_bench/create_dtoxs_splits.py --data input/dtoxs.h5ad
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

logger = logging.getLogger("create_splits")
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

    # Only assign splits for valid drug indices (small molecules)
    valid_mask = pert_idx >= 0
    valid_indices = np.where(valid_mask)[0]

    # Group valid sample indices by drug index
    drug_to_samples = {}
    for i in valid_indices:
        did = pert_idx[i]
        drug_to_samples.setdefault(did, []).append(i)

    drug_ids = sorted(drug_to_samples.keys())
    logger.info(f"Found {len(drug_ids)} unique drugs with {len(valid_indices)} total samples")
    logger.info(f"Total samples in adata: {adata.n_obs} ({np.sum(~valid_mask)} non-small-molecule excluded)")

    # Track coverage: how many folds each sample appears in test
    sample_test_count = np.zeros(adata.n_obs, dtype=int)

    for fold_i in range(1, n_folds + 1):
        col_name = f"warm-start_fold-{fold_i}"
        # Initialize all to empty string (excluded); valid drugs get train/test
        fold_values = np.full(adata.n_obs, "", dtype=object)

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

    # Verify coverage: each valid sample should be in test for >= 1 fold
    uncovered = []
    for i in valid_indices:
        if sample_test_count[i] == 0:
            uncovered.append(i)

    if uncovered:
        logger.warning(f"{len(uncovered)} samples never assigned to test! Fixing...")
        # Assign uncovered samples to test in a random fold
        for i in uncovered:
            rescue_fold = rng.integers(1, n_folds + 1)
            col_name = f"warm-start_fold-{rescue_fold}"
            adata.obs.iloc[i, adata.obs.columns.get_loc(col_name)] = "test"
            sample_test_count[i] += 1
        logger.info("Coverage fixed.")

    coverage = sample_test_count[valid_indices]
    logger.info(f"Test fold coverage: min={coverage.min()}, max={coverage.max()}, "
                f"mean={coverage.mean():.1f}, median={np.median(coverage):.0f}")

    return adata


def create_dtoxs_full_split(adata):
    """Create dtoxs_full column: all valid samples = 'train', others = ''.

    Args:
        adata: AnnData with .obs containing 'pert_idx' column

    Returns:
        adata with new .obs column 'dtoxs_full'
    """
    pert_idx = adata.obs['pert_idx'].values.astype(int)
    fold_values = np.where(pert_idx >= 0, "train", "")
    adata.obs["dtoxs_full"] = fold_values
    n_train = np.sum(fold_values == "train")
    logger.info(f"  dtoxs_full: train={n_train}, test=0, valid=0")
    return adata


def main():
    parser = argparse.ArgumentParser(description="Generate warm-start and dtoxs splits in dtoxs.h5ad")
    parser.add_argument("--data", default="03_cm_finetune_data/output/dtoxs_lincs.h5ad",
                        help="Path to the DToxS-LINCS finetune h5ad relative to PROJECT_ROOT (workflow root)")
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

    # Pre-checks
    if 'pert_idx' not in adata.obs.columns:
        logger.error("'pert_idx' column not found in adata.obs — cannot create splits")
        sys.exit(1)

    logger.info("\n=== Creating warm-start splits (10 folds, drug-stratified, no validation) ===")
    adata = create_warm_start_splits(adata, n_folds=args.n_folds, seed=args.seed)

    logger.info("\n=== Creating dtoxs full split (all-train, prediction only) ===")
    adata = create_dtoxs_full_split(adata)

    # Save
    logger.info(f"\nSaving to {data_path}...")
    adata.write_h5ad(data_path)
    logger.info("Done. New columns written:")
    for i in range(1, args.n_folds + 1):
        logger.info(f"  warm-start_fold-{i}")
    logger.info("  dtoxs_full")


if __name__ == "__main__":
    main()
