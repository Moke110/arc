#!/usr/bin/env python
"""Generate the benchmark fold columns in cm_union.h5ad (ADR 0017).

The v3 union benchmark's test sets are **byte-identical to v2's** (which are
byte-identical to v1.1's): for every fold column, the primary rows carry
v1.1's fold value verbatim (matched by the loader well id, pert_id — the same
LINCS LDS-1593 wells), and the iPSC rows are assigned under the ADR 0015
isolation rules. The two concentration-response datasets (gse244740,
gse262419) are **train-only** — they never appear in a test set:

  warm-start (10 folds, evaluation-only) — primary rows verbatim from v1.1;
    iPSC rows all-train ('train'); gse244740/gse262419 rows all-train.

  cold-drug (10 folds, evaluation-only) — primary rows verbatim from v1.1;
    iPSC rows 'train' except the fold's test drugs, excluded ('') from train;
    gse244740/gse262419 rows 'train' except the fold's test drugs, excluded
    ('') — a cold-drug test drug must not leak into train through any of the
    four datasets (gate 3 extends the v2 rule to all four).

  cm_full (1 fold, prediction-only) — every valid sample 'train'.

The benchmark therefore answers precisely: *holding the v2 test sets fixed,
does adding gse244740 + gse262419 to the training corpus change primary/iPSC
test performance?* The primary arm's test-set identity is the ADR 0017
contract: re-running the splitter on the union would drift the test sets, so
the fold columns are imported from v1.1, never regenerated.

After writing, isolation gates assert:
  1. test-set identity — for all 20 fold columns, the primary rows carry
     v1.1's fold values byte-for-byte (well-matched column equality);
  2. no shared well across train/test per fold (well identity, all four arms);
  3. cold-drug exclusion — for each cold-drug fold, the test drugs are absent
     from all four datasets' train;
  4. train-only — gse244740/gse262419 rows are 'train' (or '' under cold-drug
     exclusion) and never 'test' in any fold;
plus the ADR 0014 invariant checks over the **benchmark surface** (primary +
iPSC rows, where test rows can come from): value set subset of
{train,test,''}; warm-start drug overlap 1.0; cold-drug drug overlap 0.0;
per-fold test fraction in [0.10, 0.30]; every benchmark drug and every primary
sample in test >= 1 fold; `cm_full` = every valid sample 'train'; no stale
fold columns. Any check failing aborts with exit 1.

Usage:
    python 05_xpert_finetune_bench/dtoxs-primary/create_cm_splits.py \
        --data 03_cm_finetune_data/output/cm_union.h5ad
    python 05_xpert_finetune_bench/dtoxs-primary/create_cm_splits.py \
        --v11 ../../xpert-fda-miri-v1.1/03_cm_finetune_data/output/dtoxs_lincs.h5ad
"""

import argparse
import logging
import re
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "_lib"))
from utils import N_L1000_DRUGS, V11_DTOXS_LINCS_H5AD  # noqa: E402

logger = logging.getLogger("create_cm_splits")
logger.setLevel(logging.INFO)
if not logger.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

FULL_COL = "cm_full"
FOLD_COLS = [f"warm-start_fold-{i}" for i in range(1, 11)] + \
            [f"cold-drug_fold-{i}" for i in range(1, 11)]
SPLIT_COL_RE = re.compile(r"(fold|_full)", re.IGNORECASE)
# The benchmark surface: only these datasets can contribute test rows. The
# concentration-response datasets (gse244740, gse262419) are train-only.
BENCHMARK_DATASETS = ["dtoxs_primary", "dtoxs_ipsc"]
# The isolation rule covers every finetune dataset.
ALL_FINETUNE_DATASETS = ["dtoxs_primary", "dtoxs_ipsc",
                         "gse244740", "gse262419"]


def valid_mask(adata):
    pert_idx = adata.obs["pert_idx"].values.astype(int)
    return (pert_idx >= 0) & (pert_idx < N_L1000_DRUGS)


def is_split_column(name):
    return bool(SPLIT_COL_RE.search(name))


def drop_stale_split_columns(adata):
    stale = [c for c in adata.obs.columns if is_split_column(c)]
    if stale:
        logger.info(f"Dropping stale split columns: {sorted(stale)}")
        adata.obs.drop(columns=stale, inplace=True)
    return adata


def import_v11_folds(adata, v11_path):
    """Primary rows carry v1.1's fold values verbatim (matched by well id)."""
    if not Path(v11_path).exists():
        logger.error(f"v1.1 h5ad (fold source) not found: {v11_path}")
        sys.exit(1)
    v11 = sc.read_h5ad(v11_path)
    missing_cols = [c for c in FOLD_COLS if c not in v11.obs.columns]
    if missing_cols:
        logger.error(f"v1.1 h5ad missing fold columns: {missing_cols}")
        sys.exit(1)

    v11_well = v11.obs["pert_id"].astype(str)
    v11_map = {w: i for i, w in enumerate(v11_well)}
    prim = adata.obs["dataset_id"].astype(str) == "dtoxs_primary"
    wells = adata.obs["pert_id"].astype(str).values
    prim_rows = np.where(prim)[0]
    src_idx = np.array([v11_map.get(w, -1) for w in wells])
    if int((src_idx[prim_rows] < 0).sum()):
        n_miss = int((src_idx[prim_rows] < 0).sum())
        logger.error(f"{n_miss} primary wells not found in the v1.1 h5ad — "
                     f"cannot import fold columns")
        sys.exit(1)

    for col in FOLD_COLS:
        col_vals = np.full(adata.n_obs, "", dtype=object)
        col_vals[prim_rows] = v11.obs[col].astype(str).values[src_idx[prim_rows]]
        adata.obs[col] = col_vals
    logger.info(f"Imported {len(FOLD_COLS)} fold columns from v1.1 for "
                f"{len(prim_rows)} primary rows (matched by pert_id).")
    return adata, v11


def assign_isolation(adata):
    """iPSC + gse244740 + gse262419 rows under the train-only rules.

    warm-start: all-train for every non-primary row. cold-drug: 'train'
    except the fold's test drugs, which are excluded ('') from every dataset's
    train (gate 3 — a cold-drug test drug must not leak into train through any
    of the added datasets).
    """
    nonprim = ~(adata.obs["dataset_id"].astype(str) == "dtoxs_primary")
    nonprim_rows = np.where(nonprim)[0]
    pert_idx = adata.obs["pert_idx"].values.astype(int)

    for col in FOLD_COLS:
        family = "warm-start" if col.startswith("warm-start") else "cold-drug"
        col_vals = adata.obs[col].astype(str).values
        if family == "warm-start":
            col_vals[nonprim_rows] = "train"
        else:
            test_drugs = set(int(p) for p in pert_idx[col_vals == "test"])
            for i in nonprim_rows:
                col_vals[i] = "" if int(pert_idx[i]) in test_drugs else "train"
        adata.obs[col] = col_vals


def create_full_split(adata):
    pert_idx = adata.obs["pert_idx"].values.astype(int)
    col = np.full(adata.n_obs, "", dtype=object)
    col[valid_mask(adata)] = "train"
    adata.obs[FULL_COL] = col
    return adata


def verify_gates(adata, v11, n_folds=10):
    """The four isolation gates + the ADR 0014 invariant checks."""
    pert_idx = adata.obs["pert_idx"].values.astype(int)
    valid = valid_mask(adata)
    prim = adata.obs["dataset_id"].astype(str) == "dtoxs_primary"
    bench = adata.obs["dataset_id"].astype(str).isin(BENCHMARK_DATASETS).values
    newds = adata.obs["dataset_id"].astype(str).isin(
        ["gse244740", "gse262419"]).values
    problems = []

    # Gate 1 — test-set identity vs v1.1 (byte-level: verbatim column values
    # on the well-matched rows, not just the same well set)
    prim_wells = adata.obs.loc[prim, "pert_id"].astype(str)
    v11_wells = v11.obs["pert_id"].astype(str)
    v11_row_of = {w: i for i, w in enumerate(v11_wells)}
    prim_v11_idx = np.array([v11_row_of[w] for w in prim_wells])
    for col in FOLD_COLS:
        v3_vals = adata.obs.loc[prim, col].astype(str).values
        v11_vals = v11.obs[col].astype(str).values[prim_v11_idx]
        if not np.array_equal(v3_vals, v11_vals):
            problems.append(f"{col}: primary fold values not byte-identical "
                            f"to v1.1 (well-matched columns differ)")

    # Gate 2 — no identical well across train+test per fold (all four arms)
    for col in FOLD_COLS:
        v = adata.obs[col].astype(str).values
        train_wells = set(adata.obs["pert_id"].astype(str).values[v == "train"])
        test_wells = set(adata.obs["pert_id"].astype(str).values[v == "test"])
        shared = train_wells & test_wells
        if shared:
            problems.append(f"{col}: {len(shared)} wells in both train and test")

    # Gate 3 — cold-drug test drugs absent from all four datasets' train
    for col in FOLD_COLS:
        if not col.startswith("cold-drug"):
            continue
        v = adata.obs[col].astype(str).values
        test_drugs = set(int(p) for p in pert_idx[v == "test"])
        train_drugs = set(int(p) for p in pert_idx[v == "train"])
        leak = test_drugs & train_drugs
        if leak:
            problems.append(f"{col}: cold-drug overlap > 0.0 — "
                            f"{len(leak)} test drugs in train (any arm)")

    # Gate 4 — gse244740/gse262419 rows are train-only: never 'test'
    for col in FOLD_COLS:
        v = adata.obs[col].astype(str).values
        if int((v[newds] == "test").sum()):
            problems.append(f"{col}: {int((v[newds] == 'test').sum())} "
                            f"gse244740/gse262419 rows in a test set "
                            f"(train-only datasets must never be test)")

    # ADR 0014 invariant checks over the benchmark surface (primary + iPSC:
    # the rows a test set can come from)
    valid_drugs = sorted(set(int(p) for p in pert_idx[valid & bench]))
    target_cols = set(FOLD_COLS) | {FULL_COL}
    stale = [c for c in adata.obs.columns if is_split_column(c)
             and c not in target_cols]
    if stale:
        problems.append(f"stale split columns remain: {sorted(stale)}")

    for i in range(1, n_folds + 1):
        for col in (f"warm-start_fold-{i}", f"cold-drug_fold-{i}"):
            v = adata.obs[col].astype(str).values
            if not set(v) <= {"train", "test", ""}:
                problems.append(f"{col}: values {set(v)} outside "
                                f"{{train, test, ''}}")
            # test fraction over the benchmark surface only (the train-only
            # rows inflate n_train and would otherwise break the band)
            bv = v[bench & valid]
            n_test = int((bv == "test").sum())
            n_train = int((bv == "train").sum())
            frac = n_test / (n_test + n_train) if (n_test + n_train) else 0.0
            if not 0.10 <= frac <= 0.30:
                problems.append(f"{col}: benchmark test fraction {frac:.3f} "
                                f"outside [0.10, 0.30]")
            logger.info(f"  {col}: bench_test={n_test}, bench_train={n_train}, "
                        f"bench_test_frac={frac:.3f}")
            test_drugs = set(int(p) for p in pert_idx[v == "test"])
            train_drugs = set(int(p) for p in pert_idx[v == "train"])
            if col.startswith("warm-start"):
                missing = test_drugs - train_drugs
                if missing:
                    problems.append(f"{col}: warm-start overlap < 1.0 — "
                                    f"{len(missing)} test drugs never in train")

    # every benchmark drug in test >= 1 fold; every PRIMARY sample in test >= 1
    drug_test_count = {}
    sample_test_count = np.zeros(adata.n_obs, dtype=int)
    for col in FOLD_COLS:
        v = adata.obs[col].astype(str).values
        sample_test_count += (v == "test")
        for d in set(int(p) for p in pert_idx[v == "test"]):
            drug_test_count[d] = drug_test_count.get(d, 0) + 1
    missing_drugs = [d for d in valid_drugs if drug_test_count.get(d, 0) == 0]
    if missing_drugs:
        problems.append(f"{len(missing_drugs)} benchmark drugs never in test: "
                        f"{sorted(missing_drugs)[:10]}")
    n_uncovered_prim = int((sample_test_count[prim.values] == 0).sum())
    if n_uncovered_prim:
        problems.append(f"warm-start: {n_uncovered_prim} primary samples "
                        f"never in test")

    v = adata.obs[FULL_COL].astype(str).values
    if int((v[valid] == "train").sum()) != int(valid.sum()):
        problems.append(f"{FULL_COL}: not all valid samples are 'train'")

    if problems:
        logger.error("SPLIT GATE FAILED:")
        for p in problems:
            logger.error(f"  - {p}")
        sys.exit(1)
    logger.info("Split gates passed: primary test sets byte-identical to "
                "v1.1; no shared well across train/test; cold-drug test drugs "
                "excluded from all four datasets' train; gse244740/gse262419 "
                "train-only; warm-start overlap 1.0; benchmark test fraction "
                "in [0.10, 0.30]; every benchmark drug and every primary "
                "sample covered; no stale fold columns.")


def main():
    parser = argparse.ArgumentParser(
        description="Import v1.1 fold columns into the v3 union h5ad and "
                    "assign iPSC/gse244740/gse262419 rows under the "
                    "train-only isolation rules")
    parser.add_argument("--data", default="03_cm_finetune_data/output/cm_union.h5ad",
                        help="Union h5ad relative to the workflow root")
    parser.add_argument("--v11", default=str(V11_DTOXS_LINCS_H5AD),
                        help="v1.1 DToxS-LINCS h5ad (fold source)")
    args = parser.parse_args()

    data_path = PROJECT_ROOT / args.data
    if not data_path.exists():
        logger.error(f"Data file not found: {data_path}")
        sys.exit(1)

    logger.info(f"Loading {data_path}...")
    adata = sc.read_h5ad(data_path)
    logger.info(f"Loaded: {adata.n_obs} obs, {adata.n_vars} vars")
    if "pert_idx" not in adata.obs.columns:
        logger.error("'pert_idx' column not found in adata.obs")
        sys.exit(1)

    drop_stale_split_columns(adata)
    adata, v11 = import_v11_folds(adata, args.v11)
    assign_isolation(adata)
    create_full_split(adata)

    logger.info(f"\nSaving to {data_path}...")
    adata.write_h5ad(data_path)
    logger.info("Re-reading written h5ad for the isolation gates...")
    written = sc.read_h5ad(data_path)
    verify_gates(written, v11, n_folds=10)
    logger.info("Done. Columns written: warm-start_fold-1..10, "
                "cold-drug_fold-1..10, cm_full")


if __name__ == "__main__":
    main()
