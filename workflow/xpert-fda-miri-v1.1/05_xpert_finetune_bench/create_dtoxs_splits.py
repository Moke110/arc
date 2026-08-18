#!/usr/bin/env python
"""Generate the unified benchmark fold columns in dtoxs.h5ad.

One fold principle, applied identically to every MIRI workflow (ADR 0014).
Writes exactly these .obs columns on the finetune h5ad:

  warm-start (10 folds, evaluation-only) — drug-overlap: every drug in the
    test set also appears in the train set. Per drug, samples are randomly
    assigned ~20% (test_frac) to test and the rest to train; any sample never
    in test across the 10 folds is rescued into a random fold (keeping its
    drug >= 1 train sample there). Values: 'train' / 'test' / '' (excluded).

  cold-drug (10 folds, evaluation-only) — drug-disjoint: no drug in the test
    set appears in the train set. Per fold, drugs are selected sample-weighted
    (so the sample ratio lands near 1:4) until ~20% of all valid samples are
    test; every other drug is train. Any drug never in test is rescued into a
    random fold (capped so the fold's test fraction stays in the ratio band).
    Values: 'train' / 'test' / '' (excluded).

  dtoxs_full (1 fold, prediction-only) — every valid sample 'train'. Values:
    'train' / ''.

Ratio band: per fold, N(test) : N(train) in 1:9 - 3:7 (test_frac ~0.2).
No validation set: fold columns carry only train/test; the XPert loader's
valid -> test fallback covers the empty validation set. No category or
data-source stratification.

Only samples with a valid HG-resolved pert_idx (0 <= pert_idx < N_L1000_DRUGS)
are included in splits; non-small-molecule / excluded perturbations are left
empty (''). All prior fold/full columns are dropped first, so the h5ad ends
with exactly the 21 columns above.

After writing, a built-in invariant gate asserts and logs for every fold
column: value set subset of {train, test, ''}; warm-start drug overlap = 1.0;
cold-drug drug overlap = 0.0; per-fold test fraction in [0.10, 0.30]; every
valid drug (and every warm-start sample) in test >= 1 fold; no stale fold
columns remain. Any check failing aborts with exit 1.

Deterministic under seed=2024. Values are full strings ('train'/'test') —
never a fixed-width numpy dtype that truncates them.

Usage:
    python 05_xpert_finetune_bench/create_dtoxs_splits.py \
        --data 03_cm_finetune_data/output/dtoxs_lincs.h5ad
"""

import argparse
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "_lib"))
from utils import N_L1000_DRUGS  # noqa: E402

logger = logging.getLogger("create_dtoxs_splits")
logger.setLevel(logging.INFO)
if not logger.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

FULL_COL = "dtoxs_full"
SPLIT_COL_RE = re.compile(r"(fold|_full)", re.IGNORECASE)


def valid_mask(adata):
    pert_idx = adata.obs["pert_idx"].values.astype(int)
    return (pert_idx >= 0) & (pert_idx < N_L1000_DRUGS)


def _group_samples(adata):
    pert_idx = adata.obs["pert_idx"].values.astype(int)
    valid = valid_mask(adata)
    valid_idx = np.where(valid)[0]
    drug_to_samples = defaultdict(list)
    for i in valid_idx:
        drug_to_samples[int(pert_idx[i])].append(int(i))
    drug_ids = sorted(drug_to_samples.keys())
    counts = {d: len(drug_to_samples[d]) for d in drug_ids}
    return pert_idx, valid, valid_idx, drug_to_samples, drug_ids, counts


def create_warm_start_splits(adata, n_folds=10, test_frac=0.2, seed=2024):
    """Create warm-start (drug-overlap) folds: ~test_frac of each drug's
    samples to test per fold, the rest train; sample coverage rescued."""
    rng = np.random.default_rng(seed)
    pert_idx, _, valid_idx, drug_to_samples, drug_ids, _ = _group_samples(adata)
    sample_test_count = np.zeros(adata.n_obs, dtype=int)
    fold_cols = {}
    for fold_i in range(1, n_folds + 1):
        col = np.full(adata.n_obs, "", dtype=object)
        for did in drug_ids:
            samples = drug_to_samples[did]
            n_test = max(1, int(np.round(len(samples) * test_frac)))
            shuffled = rng.permutation(samples)
            test_idx = shuffled[:n_test]
            train_idx = shuffled[n_test:]
            col[test_idx] = "test"
            col[train_idx] = "train"
            sample_test_count[test_idx] += 1
        fold_cols[f"warm-start_fold-{fold_i}"] = col
    uncovered = valid_idx[sample_test_count[valid_idx] == 0]
    for i in uncovered:
        did = int(pert_idx[i])
        samples = drug_to_samples[did]
        for _ in range(n_folds):
            rescue_fold = int(rng.integers(1, n_folds + 1))
            col = fold_cols[f"warm-start_fold-{rescue_fold}"]
            # Flipping this sample to test must leave its drug >= 1 train
            # sample here, else the fold's drug-overlap invariant breaks.
            # Only n < 3-sample drugs (which the ratio band cannot represent)
            # can violate it; they fall through and the gate reports them.
            if int(np.sum(col[samples] == "train")) >= 2:
                col[i] = "test"
                break
    for name, col in fold_cols.items():
        adata.obs[name] = col
    return adata


def create_cold_drug_splits(adata, n_folds=10, test_frac=0.2, seed=2024):
    """Create cold-drug (drug-disjoint) folds: sample-weighted drug selection
    to ~test_frac of all valid samples as test per fold; drug coverage
    rescued."""
    rng = np.random.default_rng(seed)
    _, _, _, drug_to_samples, drug_ids, counts = _group_samples(adata)
    total = float(sum(counts.values()))
    target = test_frac * total

    test_drugs_per_fold = []
    for _ in range(n_folds):
        remaining = list(drug_ids)
        remaining_counts = [counts[d] for d in remaining]
        test_drugs = set()
        test_count = 0.0
        while test_count < target and remaining:
            w = np.asarray(remaining_counts, dtype=float)
            w /= w.sum()
            choice = int(rng.choice(len(remaining), p=w))
            d = remaining.pop(choice)
            remaining_counts.pop(choice)
            test_drugs.add(d)
            test_count += counts[d]
        test_drugs_per_fold.append(test_drugs)

    def fold_frac(fold_set):
        return sum(counts[x] for x in fold_set) / total

    covered = defaultdict(int)
    for s in test_drugs_per_fold:
        for d in s:
            covered[d] += 1
    uncovered_drugs = [d for d in drug_ids if covered[d] == 0]
    rng.shuffle(uncovered_drugs)
    for d in uncovered_drugs:
        candidates = list(range(n_folds))
        rng.shuffle(candidates)
        chosen = None
        for f in candidates:
            if fold_frac(test_drugs_per_fold[f] | {d}) <= 0.30:
                chosen = f
                break
        if chosen is None:
            chosen = min(range(n_folds), key=lambda f: fold_frac(test_drugs_per_fold[f]))
        test_drugs_per_fold[chosen].add(d)
        covered[d] += 1

    for fold_i, test_drugs in enumerate(test_drugs_per_fold, start=1):
        col = np.full(adata.n_obs, "", dtype=object)
        for d in test_drugs:
            col[drug_to_samples[d]] = "test"
        for d in drug_ids:
            if d not in test_drugs:
                col[drug_to_samples[d]] = "train"
        adata.obs[f"cold-drug_fold-{fold_i}"] = col
    return adata


def create_full_split(adata):
    """Create the all-train (prediction-only) column: every valid sample 'train'."""
    pert_idx = adata.obs["pert_idx"].values.astype(int)
    col = np.full(adata.n_obs, "", dtype=object)
    col[pert_idx >= 0] = "train"
    adata.obs[FULL_COL] = col
    return adata


def is_split_column(name):
    return bool(SPLIT_COL_RE.search(name))


def drop_stale_split_columns(adata):
    stale = [c for c in adata.obs.columns if is_split_column(c)]
    if stale:
        logger.info(f"Dropping stale split columns: {sorted(stale)}")
        adata.obs.drop(columns=stale, inplace=True)
    return adata


def verify_split_gate(adata, n_folds=10):
    """Invariant gate on the written fold columns. Any failure exits 1."""
    pert_idx, _, valid_idx, _, _, _ = _group_samples(adata)
    valid_drugs = sorted(set(int(p) for p in pert_idx[valid_idx]))
    problems = []

    target_cols = {f"warm-start_fold-{i}" for i in range(1, n_folds + 1)} \
        | {f"cold-drug_fold-{i}" for i in range(1, n_folds + 1)} \
        | {FULL_COL}
    stale = [c for c in adata.obs.columns if is_split_column(c) and c not in target_cols]
    if stale:
        problems.append(f"stale split columns remain: {sorted(stale)}")

    for i in range(1, n_folds + 1):
        for col in (f"warm-start_fold-{i}", f"cold-drug_fold-{i}"):
            v = adata.obs[col].astype(str).values
            if not set(v) <= {"train", "test", ""}:
                problems.append(f"{col}: values {set(v)} outside "
                                f"{{train, test, ''}}")
            test_mask = v == "test"
            train_mask = v == "train"
            n_test = int(test_mask.sum())
            n_train = int(train_mask.sum())
            frac = n_test / (n_test + n_train) if (n_test + n_train) else 0.0
            if not 0.10 <= frac <= 0.30:
                problems.append(f"{col}: test fraction {frac:.3f} outside "
                                f"[0.10, 0.30] (n_test={n_test}, n_train={n_train})")
            logger.info(f"  {col}: test={n_test}, train={n_train}, "
                        f"test_frac={frac:.3f}")
            test_drugs = set(int(p) for p in pert_idx[test_mask])
            train_drugs = set(int(p) for p in pert_idx[train_mask])
            if col.startswith("warm-start"):
                missing = test_drugs - train_drugs
                if missing:
                    problems.append(f"{col}: warm-start drug overlap < 1.0 — "
                                    f"{len(missing)} test drugs never in train")
            else:
                inter = test_drugs & train_drugs
                if inter:
                    problems.append(f"{col}: cold-drug overlap > 0.0 — "
                                    f"{len(inter)} drugs in both test and train")

    drug_test_count = defaultdict(int)
    for i in range(1, n_folds + 1):
        for family in ("warm-start", "cold-drug"):
            v = adata.obs[f"{family}_fold-{i}"].astype(str).values
            for d in set(int(p) for p in pert_idx[v == "test"]):
                drug_test_count[(family, d)] += 1
    for family in ("warm-start", "cold-drug"):
        missing = [d for d in valid_drugs if drug_test_count[(family, d)] == 0]
        if missing:
            problems.append(f"{family}: {len(missing)} drugs never in test: "
                            f"{sorted(missing)[:10]}")

    sample_test_count = np.zeros(adata.n_obs, dtype=int)
    for i in range(1, n_folds + 1):
        sample_test_count += adata.obs[f"warm-start_fold-{i}"].astype(str).values == "test"
    n_uncovered = int((sample_test_count[valid_idx] == 0).sum())
    if n_uncovered:
        problems.append(f"warm-start: {n_uncovered} valid samples never in test")

    v = adata.obs[FULL_COL].astype(str).values
    if int((v[valid_idx] == "train").sum()) != len(valid_idx):
        problems.append(f"{FULL_COL}: not all valid samples are 'train'")

    if problems:
        logger.error("SPLIT GATE FAILED:")
        for p in problems:
            logger.error(f"  - {p}")
        sys.exit(1)
    logger.info("Split gate passed: values in {train, test, ''}; warm-start "
                "drug overlap 1.0; cold-drug drug overlap 0.0; per-fold test "
                "fraction in [0.10, 0.30]; every drug/sample covered; no stale "
                "fold columns.")


def main():
    parser = argparse.ArgumentParser(
        description="Generate the unified warm-start + cold-drug + full split "
                    "columns in dtoxs.h5ad")
    parser.add_argument("--data", default="03_cm_finetune_data/output/dtoxs_lincs.h5ad",
                        help="Path to dtoxs.h5ad relative to the workflow root")
    parser.add_argument("--seed", type=int, default=2024, help="Random seed")
    parser.add_argument("--n_folds", type=int, default=10,
                        help="Number of warm-start / cold-drug folds")
    args = parser.parse_args()

    data_path = PROJECT_ROOT / args.data
    if not data_path.exists():
        logger.error(f"Data file not found: {data_path}")
        sys.exit(1)

    logger.info(f"Loading {data_path}...")
    adata = sc.read_h5ad(data_path)
    logger.info(f"Loaded: {adata.n_obs} obs, {adata.n_vars} vars")

    if "pert_idx" not in adata.obs.columns:
        logger.error("'pert_idx' column not found in adata.obs — cannot "
                     "create splits")
        sys.exit(1)

    drop_stale_split_columns(adata)

    logger.info("\n=== Creating warm-start splits "
                f"({args.n_folds} folds, drug-overlap, no validation) ===")
    adata = create_warm_start_splits(adata, n_folds=args.n_folds,
                                     test_frac=0.2, seed=args.seed)

    logger.info("\n=== Creating cold-drug splits "
                f"({args.n_folds} folds, drug-disjoint, no validation) ===")
    adata = create_cold_drug_splits(adata, n_folds=args.n_folds,
                                    test_frac=0.2, seed=args.seed)

    logger.info("\n=== Creating full split (all-train, prediction only) ===")
    adata = create_full_split(adata)

    logger.info(f"\nSaving to {data_path}...")
    adata.write_h5ad(data_path)
    logger.info("Re-reading written h5ad for the invariant gate...")
    written = sc.read_h5ad(data_path)
    verify_split_gate(written, n_folds=args.n_folds)
    logger.info("Done. Columns written: "
                f"warm-start_fold-1..{args.n_folds}, "
                f"cold-drug_fold-1..{args.n_folds}, {FULL_COL}")


if __name__ == "__main__":
    main()
