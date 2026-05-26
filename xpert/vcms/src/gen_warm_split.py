"""
Generate warm 7-fold splits for dtoxs dataset.

- train:valid:test = 7:1.5:1.5
- Each sample must appear in train at least once AND in test at least once across all 7 folds.
- Writes warm_fold-{1..7} columns to dtoxs.h5ad and updates dtoxs.yaml.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata
import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

H5AD_PATH = PROJECT_ROOT / "finetune" / "data" / "dtoxs.h5ad"
YAML_PATH = PROJECT_ROOT / "finetune" / "config" / "dtoxs.yaml"

N_FOLDS = 7
TRAIN_RATIO = 0.70
VALID_RATIO = 0.15
TEST_RATIO = 0.15
SEED = 42


def main():
    rng = np.random.default_rng(SEED)

    # --- Load data ---
    print(f"Loading {H5AD_PATH}...")
    adata = anndata.read_h5ad(H5AD_PATH)
    n_samples = adata.n_obs
    print(f"  {n_samples} samples")

    # --- Remove old split columns ---
    old_cols = [c for c in adata.obs.columns if any(
        prefix in c for prefix in ("split_", "random_fold-", "cold-drug_fold-", "cold-cate_fold-")
    )]
    print(f"Removing {len(old_cols)} old split columns: {old_cols}")
    adata.obs.drop(columns=old_cols, inplace=True)

    # --- Generate initial random splits ---
    print(f"Generating {N_FOLDS} folds ({TRAIN_RATIO}:{VALID_RATIO}:{TEST_RATIO})...")
    indices = np.arange(n_samples)
    fold_assignments = np.zeros((N_FOLDS, n_samples), dtype=np.int8)  # 0=train, 1=valid, 2=test

    for k in range(N_FOLDS):
        shuffled = rng.permutation(indices)
        n_train = int(n_samples * TRAIN_RATIO)
        n_valid = int(n_samples * VALID_RATIO)
        # n_test = n_samples - n_train - n_valid

        train_idx = shuffled[:n_train]
        valid_idx = shuffled[n_train : n_train + n_valid]
        test_idx = shuffled[n_train + n_valid :]

        fold_assignments[k, train_idx] = 0
        fold_assignments[k, valid_idx] = 1
        fold_assignments[k, test_idx] = 2

        print(f"  Fold {k+1}: train={len(train_idx)} valid={len(valid_idx)} test={len(test_idx)}")

    # --- Coverage check and fix ---
    in_train = np.sum(fold_assignments == 0, axis=0)  # per sample
    in_test = np.sum(fold_assignments == 2, axis=0)

    print(f"\nInitial coverage:")
    print(f"  in_train=0: {np.sum(in_train == 0)} samples")
    print(f"  in_test=0:  {np.sum(in_test == 0)} samples")

    # Fix in_train == 0
    zero_train = np.where(in_train == 0)[0]
    if len(zero_train) > 0:
        print(f"\nFixing {len(zero_train)} samples with in_train=0...")
        for sample_i in zero_train:
            # Find a fold where this sample is NOT in train (it must be, since in_train=0)
            # Pick a fold at random, swap with a train sample that has in_train > 1
            k = rng.integers(0, N_FOLDS)
            fold_k = fold_assignments[k, :]
            train_in_fold = np.where(fold_k == 0)[0]
            # Among train samples in this fold, find those with in_train > 1
            candidates = train_in_fold[in_train[train_in_fold] > 1]
            if len(candidates) == 0:
                # Try another fold
                for alt_k in range(N_FOLDS):
                    fold_alt = fold_assignments[alt_k, :]
                    train_in_alt = np.where(fold_alt == 0)[0]
                    candidates = train_in_alt[in_train[train_in_alt] > 1]
                    if len(candidates) > 0:
                        k = alt_k
                        break
            if len(candidates) == 0:
                print(f"  WARNING: could not fix sample {sample_i}")
                continue
            swap_j = rng.choice(candidates)
            # Swap assignments in this fold
            fold_assignments[k, sample_i], fold_assignments[k, swap_j] = (
                fold_assignments[k, swap_j], fold_assignments[k, sample_i],
            )
            # Update counts
            in_train[sample_i] += 1
            in_train[swap_j] -= 1

    # Fix in_test == 0
    zero_test = np.where(in_test == 0)[0]
    if len(zero_test) > 0:
        print(f"Fixing {len(zero_test)} samples with in_test=0...")
        for sample_i in zero_test:
            k = rng.integers(0, N_FOLDS)
            fold_k = fold_assignments[k, :]
            test_in_fold = np.where(fold_k == 2)[0]
            candidates = test_in_fold[in_test[test_in_fold] > 1]
            if len(candidates) == 0:
                for alt_k in range(N_FOLDS):
                    fold_alt = fold_assignments[alt_k, :]
                    test_in_alt = np.where(fold_alt == 2)[0]
                    candidates = test_in_alt[in_test[test_in_alt] > 1]
                    if len(candidates) > 0:
                        k = alt_k
                        break
            if len(candidates) == 0:
                print(f"  WARNING: could not fix sample {sample_i}")
                continue
            swap_j = rng.choice(candidates)
            fold_assignments[k, sample_i], fold_assignments[k, swap_j] = (
                fold_assignments[k, swap_j], fold_assignments[k, sample_i],
            )
            in_test[sample_i] += 1
            in_test[swap_j] -= 1

    # --- Final check ---
    in_train = np.sum(fold_assignments == 0, axis=0)
    in_test = np.sum(fold_assignments == 2, axis=0)
    print(f"\nFinal coverage:")
    print(f"  in_train=0: {np.sum(in_train == 0)} samples")
    print(f"  in_test=0:  {np.sum(in_test == 0)} samples")
    print(f"  in_train: min={in_train.min()} max={in_train.max()} mean={in_train.mean():.1f}")
    print(f"  in_test:  min={in_test.min()} max={in_test.max()} mean={in_test.mean():.1f}")

    # --- Write to h5ad ---
    label_map = {0: "train", 1: "valid", 2: "test"}
    for k in range(N_FOLDS):
        col_name = f"warm_fold-{k+1}"
        labels = [label_map[v] for v in fold_assignments[k, :]]
        adata.obs[col_name] = pd.Categorical(labels, categories=["train", "valid", "test"])

    print(f"\nWriting {H5AD_PATH}...")
    adata.write_h5ad(H5AD_PATH)
    print("  Done.")

    # --- Update dtoxs.yaml ---
    print(f"\nUpdating {YAML_PATH}...")
    with open(YAML_PATH) as f:
        config = yaml.safe_load(f)

    config["split_methods"] = {
        "warm": {"folds": [f"warm_fold-{k+1}" for k in range(N_FOLDS)]},
    }

    with open(YAML_PATH, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print("  Done.")
    print(f"\nNew split_methods: {config['split_methods']}")


if __name__ == "__main__":
    main()
