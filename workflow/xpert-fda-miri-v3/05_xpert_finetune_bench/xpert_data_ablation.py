#!/usr/bin/env python3
"""Step 05 — Data-ablation benchmark for the v3 MIRI finetune (ADR 0018).

Measures how the finetune-corpus SCALE affects XPert test performance in
warm-start and cold-drug modes, **decoupled from the official dtoxs-primary
folds**. The experiment is a diagnostic analysis only: it changes no model
used downstream (step 06 cm_full and steps 07-11 are untouched).

Design (3 seeds x 2 modes x 5 fractions = 30 finetune runs + 3 anchors):
  - Each seed holds one **unified test set** (~15% of the 4-dataset union,
    uniform random sample-wise, composition mirroring the union), shared by
    every fraction of both modes.
  - Train pools: warm-start = union minus test; cold-drug = union minus test
    minus every sample of the test drugs across all four datasets (strictly
    drug-disjoint at every fraction).
  - Fractions {1/16, 1/8, 1/4, 1/2, 1.0}: each a dataset-stratified
    f-subsample of the seed's pool (per-dataset row counts x f, composition
    held constant), drawn with a derived deterministic seed per
    (seed, mode, fraction). 1.0 is re-run on this scheme (never the
    dtoxs-primary bench).
  - Training recipe byte-identical to the official fold bench (4 epochs,
    lr 0.002, batch 32, cm.yaml, unimol+HG drug features, step-04 L1000
    pretrain init). The L1000-pretrain weight anchors the curve at 0 rows
    (3 zero-shot evaluations, no training).

Run modes:
  --dry-run   build the masks, run the isolation gates, persist the mask
              manifest, and exit — no GPU work.
  (default)   sequential 30-run training + per-run test evaluation (strictly
              one training process at a time), then the 3 pretrain anchors,
              then aggregate -> plot -> report. Idempotent per
              (seed, mode, fraction) via caches; --force redoes finished
              runs; an interrupted job resumes from where it stopped.

Outputs under 05_xpert_finetune_bench/data-ablation/:
  ablation_masks.csv   per-(seed, mode, fraction) test/train well sets
                       (+ n_train/n_test/n_test_drugs/n_absent_drugs)
  work/                temporary h5ads (union copy + ablation fold columns);
                       cm_union.h5ad and its 21 official fold columns are
                       never written.
  weights/             per-run best weights
  results/             per-run test result CSVs
  ablation.csv         33 rows (30 runs + 3 pretrain anchors)
  ablation.png         2 rows x 4 cols (rows warm-start | cold-drug;
                       columns RMSE | PCC | Spearman | R2; x = fraction on a
                       log scale; 3 seed points + mean +/- std; pretrain
                       anchor at fraction -> 0)
  ablation_report.md   scale-curve + per-fraction drug-coverage summary

Usage:
    python 05_xpert_finetune_bench/xpert_data_ablation.py --dry-run
    python 05_xpert_finetune_bench/xpert_data_ablation.py
"""

import argparse
import copy
import hashlib
import shutil
import sys
import yaml
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT  # noqa: E402
from xpert_finetune_bench import (  # noqa: E402
    is_test_csv,
    run_xpert,
    STAGING_ROOT,
    DEFAULT_PRETRAINED,
)

STEP_DIR = WORKFLOW_ROOT / "05_xpert_finetune_bench"
OUTPUT_DIR = STEP_DIR / "data-ablation"
WORK_DIR = OUTPUT_DIR / "work"
WEIGHTS_DIR = OUTPUT_DIR / "weights"
RESULTS_DIR = OUTPUT_DIR / "results"
CONFIG_DIR = STEP_DIR / "config"
UNION_H5AD = WORKFLOW_ROOT / "03_cm_finetune_data" / "output" / "cm_union.h5ad"
FINETUNE_CONFIG = WORKFLOW_ROOT / "06_xpert_finetune" / "config" / "cm.yaml"
PRETRAIN_EVAL_CONFIG = STEP_DIR / "config" / "cm_pretrain_eval.yaml"

SEEDS = [1, 2, 3]
MODES = ["warm-start", "cold-drug"]
FRACTIONS = [1 / 16, 1 / 8, 1 / 4, 1 / 2, 1.0]
FRAC_TAGS = ["1of16", "1of8", "1of4", "1of2", "full"]
TEST_FRACTION = 0.15
TEST_FRACTION_TOL = 0.01
NUM_EPOCHS = 4
MASK_CSV = OUTPUT_DIR / "ablation_masks.csv"
METRICS = ["MSE_deg", "RMSE_deg", "MAE_deg", "Pearson_deg", "Spearman_deg",
           "R2_deg"]
# (panel label, metric column) — the four headline panels per mode.
PLOT_METRICS = [("RMSE", "RMSE_deg"), ("PCC", "Pearson_deg"),
                ("Spearman", "Spearman_deg"), ("R2", "R2_deg")]
ANCHOR_LABEL = "anchor"
# Canonical column order of the mask manifest (ablation_masks.csv).
MASK_COLUMNS = ["seed", "mode", "fraction", "fraction_tag", "n_test",
                "n_train", "n_test_drugs", "n_absent_drugs", "test_wells",
                "train_wells"]


def encode_wells(wells):
    return ";".join(wells)


def decode_wells(s):
    return set(s.split(";"))


def log(msg, level="INFO"):
    print(f"[{level}] {msg}")


def die(msg):
    log(msg, "ERROR")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Split construction
# ---------------------------------------------------------------------------
def derived_seed(seed, mode, fraction):
    """Deterministic derived RNG seed per (seed, mode, fraction)."""
    h = hashlib.md5(f"{seed}|{mode}|{fraction}".encode()).hexdigest()
    return int(h, 16)


def build_test_indices(adata, seed):
    """One unified test set (~15% of the union, uniform random per seed).

    The draw rescues any drug whose ALL samples landed in test by moving one
    sample back to the pool, guaranteeing warm-start f=1 train-pool overlap 1.0.
    """
    rng = np.random.default_rng(seed)
    n = adata.n_obs
    n_test = int(round(TEST_FRACTION * n))
    test = set(rng.choice(n, size=n_test, replace=False).tolist())
    pert = adata.obs["pert_idx"].astype(int).values
    by_drug = defaultdict(list)
    for i in range(n):
        by_drug[pert[i]].append(i)
    for d in sorted({pert[i] for i in test}):
        samples = by_drug[d]
        if all(i in test for i in samples):
            test.remove(samples[-1])
    return test


def _stratify_sub(pool, ds, fraction, rng):
    """Dataset-stratified f-subsample of a pool (composition held constant).

    Per-dataset counts are allocated by largest-remainder (Hamilton)
    apportionment so each dataset's share of the subsample deviates from its
    pool share by at most 1/n_train, and the total is exactly
    round(N * fraction). Deterministic: datasets in sorted order, ties in the
    remainder allocation broken by dataset name.
    """
    ds_pool = ds[pool]
    names = np.unique(ds_pool)  # sorted, deterministic
    counts = {d: int((ds_pool == d).sum()) for d in names}
    n = len(pool)
    k_total = int(round(n * fraction))
    if k_total <= 0:
        return frozenset()
    exact = {d: counts[d] * k_total / n for d in names}
    seats = {d: int(exact[d]) for d in names}
    rem = {d: exact[d] - seats[d] for d in names}
    shortfall = k_total - sum(seats.values())
    for d in sorted(names, key=lambda d: (-rem[d], d))[:shortfall]:
        seats[d] += 1
    chosen = []
    for d in names:
        k = seats[d]
        if k > 0:
            idxs = pool[ds_pool == d]
            chosen.extend(rng.choice(idxs, size=k, replace=False).tolist())
    return frozenset(chosen)


def build_all_masks(adata):
    """Per-(seed, mode, fraction) masks for all seeds. Deterministic."""
    pert = adata.obs["pert_idx"].astype(int).values
    ds = adata.obs["dataset_id"].astype(str).values
    wells = adata.obs["well"].astype(str).values
    n = adata.n_obs
    idx = np.arange(n)

    rows = []
    for seed in SEEDS:
        test = build_test_indices(adata, seed)
        in_test = np.zeros(n, dtype=bool)
        in_test[list(test)] = True
        test_drugs = set(int(p) for p in pert[in_test])
        warm_pool = idx[~in_test]
        keep = ~np.isin(pert, np.array(sorted(test_drugs)))
        cold_pool = idx[keep]
        for mode in MODES:
            pool = warm_pool if mode == "warm-start" else cold_pool
            for frac, tag in zip(FRACTIONS, FRAC_TAGS):
                rng = np.random.default_rng(derived_seed(seed, mode, frac))
                train = _stratify_sub(pool, ds, frac, rng)
                train_drugs = set(int(p) for p in pert[sorted(train)])
                rows.append({
                    "seed": seed,
                    "mode": mode,
                    "fraction": frac,
                    "fraction_tag": tag,
                    "n_test": len(test),
                    "n_train": len(train),
                    "n_test_drugs": len(test_drugs),
                    "n_absent_drugs": len(test_drugs - train_drugs),
                    "test_wells": encode_wells(wells[sorted(test)]),
                    "train_wells": encode_wells(wells[sorted(train)]),
                })
    df = pd.DataFrame(rows, columns=MASK_COLUMNS)
    return df


# ---------------------------------------------------------------------------
# Isolation gates
# ---------------------------------------------------------------------------
def _well_sets(mask_row):
    return decode_wells(mask_row["test_wells"]), decode_wells(
        mask_row["train_wells"])


def _share_of(ds, idxs):
    if len(idxs) == 0:
        return {}
    c = pd.Series(ds[idxs]).value_counts(normalize=True)
    return {str(d): float(v) for d, v in c.items()}


def verify_gates(adata, masks_df):
    """Assert the ADR 0018 design invariants; exit 1 on any failure."""
    pert = adata.obs["pert_idx"].astype(int).values
    ds = adata.obs["dataset_id"].astype(str).values
    well_to_idx = {w: i for i, w in
                   enumerate(adata.obs["well"].astype(str).values)}
    n = adata.n_obs
    problems = []

    for seed in SEEDS:
        sub = masks_df[masks_df["seed"] == seed]
        test_wells, _ = _well_sets(sub.iloc[0])
        test_frac = sub.iloc[0]["n_test"] / n
        if abs(test_frac - TEST_FRACTION) > TEST_FRACTION_TOL:
            problems.append(f"seed {seed}: test fraction {test_frac:.4f} not "
                            f"≈{TEST_FRACTION}")
        test_well_idx = {well_to_idx[w] for w in test_wells}
        in_test_mask = np.zeros(n, dtype=bool)
        in_test_mask[sorted(test_well_idx)] = True
        test_drugs = set(int(p) for p in pert[in_test_mask])

        # warm-start f=1 train pool must contain every test drug (overlap 1.0)
        ws_full = sub[(sub["mode"] == "warm-start") &
                      (sub["fraction"] == 1.0)]
        if len(ws_full) != 1 or ws_full.iloc[0]["n_absent_drugs"] != 0:
            problems.append(f"seed {seed}: warm-start f=1 overlap < 1.0 — "
                            f"{ws_full.iloc[0]['n_absent_drugs'] if len(ws_full) else '?'} "
                            f"test drugs absent from the train pool")

        for _, row in sub.iterrows():
            tw, tr = _well_sets(row)
            shared = tw & tr
            if shared:
                problems.append(f"seed {row['seed']} {row['mode']} "
                                f"{row['fraction_tag']}: {len(shared)} wells "
                                f"in both test and train")
            if row["n_train"] != len(tr) or row["n_test"] != len(tw):
                problems.append(f"seed {row['seed']} {row['mode']} "
                                f"{row['fraction_tag']}: n_train/n_test "
                                f"({row['n_train']}/{row['n_test']}) mismatch "
                                f"the persisted well sets "
                                f"({len(tr)}/{len(tw)})")

            # cold-drug: test drugs absent from train at EVERY fraction
            if row["mode"] == "cold-drug":
                if row["n_absent_drugs"] != row["n_test_drugs"]:
                    problems.append(f"seed {seed} cold-drug "
                                    f"{row['fraction_tag']}: overlap > 0.0 — "
                                    f"{row['n_test_drugs'] - row['n_absent_drugs']} "
                                    f"test drugs leaked into train")

            # dataset-stratified composition stability vs the mode's pool.
            # Apportionment bounds each dataset's share drift by 1/n_train, so
            # the tolerance scales with the subsample size (tightens at scale).
            if row["mode"] == "warm-start":
                pool = np.where(~in_test_mask)[0]
            else:
                keep = ~np.isin(pert, np.array(sorted(test_drugs)))
                pool = np.where(keep)[0]
            pool_share = _share_of(ds, pool)
            tr_idx = np.array(sorted(well_to_idx[w] for w in tr))
            sub_share = _share_of(ds, tr_idx)
            comp_tol = 1.0 / max(1, row["n_train"]) + 0.01
            for d in set(pool_share) | set(sub_share):
                delta = abs(pool_share.get(d, 0.0) - sub_share.get(d, 0.0))
                if delta > comp_tol:
                    problems.append(f"seed {seed} {row['mode']} "
                                    f"{row['fraction_tag']}: dataset {d} "
                                    f"share {sub_share.get(d, 0.0):.3f} "
                                    f"drifted from pool "
                                    f"{pool_share.get(d, 0.0):.3f} "
                                    f"(tol {comp_tol:.3f})")

        # warm-start f=1 train set must be the full pool (1.0 re-run scheme)
        ws_full = sub[(sub["mode"] == "warm-start") &
                      (sub["fraction"] == 1.0)].iloc[0]
        _, ws_tr = _well_sets(ws_full)
        if len(ws_tr) != n - ws_full["n_test"]:
            problems.append(f"seed {seed}: warm-start f=1 train size "
                            f"{len(ws_tr)} != union minus test "
                            f"{n - ws_full['n_test']}")

    if problems:
        log("SPLIT GATE FAILED:", "ERROR")
        for p in problems:
            log(f"  - {p}", "ERROR")
        sys.exit(1)
    log(f"Split gates passed for {len(SEEDS)} seeds: test size ≈{TEST_FRACTION}, "
        f"train/test well-disjoint, warm-start f=1 overlap 1.0, cold-drug "
        f"overlap 0.0 at every fraction, dataset-stratified composition "
        f"stable within the apportionment bound (1/n_train).")
    return True


# ---------------------------------------------------------------------------
# Mask manifest persistence + byte-stability
# ---------------------------------------------------------------------------
def persist_masks(fresh, force):
    if MASK_CSV.exists() and not force:
        existing = pd.read_csv(MASK_CSV)
        if masks_equal(fresh, existing):
            log(f"Mask manifest present and byte-stable: {MASK_CSV} "
                f"({len(existing)} rows)")
            return existing
        die("Mask manifest on disk differs from the re-derived masks for the "
            "same seeds — refuse to overwrite silently. Pass --force to "
            "regenerate.")
    fresh.to_csv(MASK_CSV, index=False)
    log(f"Wrote mask manifest: {MASK_CSV} ({len(fresh)} rows)")
    return fresh


def masks_equal(a, b):
    if list(a.columns) != MASK_COLUMNS or list(b.columns) != MASK_COLUMNS:
        return False
    a, b = a.copy(), b.copy()
    a["fraction"] = a["fraction"].astype(float).round(8)
    b["fraction"] = b["fraction"].astype(float).round(8)
    a["seed"] = a["seed"].astype(int)
    b["seed"] = b["seed"].astype(int)
    a["n_test"] = a["n_test"].astype(int)
    b["n_test"] = b["n_test"].astype(int)
    a["n_train"] = a["n_train"].astype(int)
    b["n_train"] = b["n_train"].astype(int)
    a["n_test_drugs"] = a["n_test_drugs"].astype(int)
    b["n_test_drugs"] = b["n_test_drugs"].astype(int)
    a["n_absent_drugs"] = a["n_absent_drugs"].astype(int)
    b["n_absent_drugs"] = b["n_absent_drugs"].astype(int)
    return a.equals(b)


# ---------------------------------------------------------------------------
# Temporary h5ads + run configs (never touch cm_union.h5ad)
# ---------------------------------------------------------------------------
def prepare_staging(adata, masks_df):
    well_idx = {w: i for i, w in enumerate(adata.obs["well"].astype(str).values)}
    n = adata.n_obs
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)

    for seed in SEEDS:
        for mode in MODES:
            sub = masks_df[(masks_df["seed"] == seed) &
                           (masks_df["mode"] == mode)]
            out = adata.copy()
            for _, row in sub.iterrows():
                col = np.full(n, "", dtype=object)
                for w in decode_wells(row["train_wells"]):
                    col[well_idx[w]] = "train"
                for w in decode_wells(row["test_wells"]):
                    col[well_idx[w]] = "test"
                out.obs[f"ablation_{row['fraction_tag']}"] = col
            path = WORK_DIR / f"seed{seed}_{mode}.h5ad"
            out.write_h5ad(path)
            log(f"Wrote temp h5ad: {path} ({out.n_obs} obs)")

        # pretrain anchor temp h5ad: seed test wells marked test
        sub = masks_df[masks_df["seed"] == seed]
        out = adata.copy()
        col = np.full(n, "", dtype=object)
        for w in decode_wells(sub.iloc[0]["test_wells"]):
            col[well_idx[w]] = "test"
        out.obs["ablation_test"] = col
        path = WORK_DIR / f"seed{seed}_pretrain.h5ad"
        out.write_h5ad(path)
        log(f"Wrote pretrain temp h5ad: {path}")

    base = yaml.safe_load(Path(FINETUNE_CONFIG).read_text())
    for seed in SEEDS:
        for mode in MODES:
            for frac, tag in zip(FRACTIONS, FRAC_TAGS):
                cfg = copy.deepcopy(base)
                cfg["dataset"]["DTOX_data_root"] = (
                    f"05_xpert_finetune_bench/data-ablation/work/"
                    f"seed{seed}_{mode}.h5ad")
                key = f"ablation_seed{seed}_{mode}_{tag}"
                cfg["split_methods"] = {key: {"folds": [f"ablation_{tag}"]}}
                (CONFIG_DIR / f"{key}.yaml").write_text(yaml.safe_dump(cfg))
    base_pe = yaml.safe_load(Path(PRETRAIN_EVAL_CONFIG).read_text())
    for seed in SEEDS:
        cfg = copy.deepcopy(base_pe)
        cfg["dataset"]["DTOX_data_root"] = (
            f"05_xpert_finetune_bench/data-ablation/work/"
            f"seed{seed}_pretrain.h5ad")
        key = f"ablation_pretrain_seed{seed}"
        cfg["split_methods"] = {key: {"folds": ["ablation_test"]}}
        (CONFIG_DIR / f"{key}.yaml").write_text(yaml.safe_dump(cfg))


def _split_key(seed, mode, tag):
    return f"ablation_seed{seed}_{mode}_{tag}"


def _model_name(split_key):
    # The fold-1 suffix mirrors src.train's fixed output naming
    # ({dataset}_{split_method}_fold-{fold_idx}), which the test-mode weight
    # lookup depends on; it is not a real fold.
    return f"cm_{split_key}_fold-1"


def _eval_and_collect(key, staging, device, config, extra):
    """Run src.train --mode test, then move the per-fold result CSV to
    RESULTS_DIR. Shared by the finetune runs and the pretrain anchors."""
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True, exist_ok=True)
    run_xpert("test", key, device, config=config, extra=extra)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    moved = False
    for csv in staging.glob("*_train_result.csv"):
        shutil.move(str(csv), str(RESULTS_DIR / csv.name))
        moved = True
    shutil.rmtree(staging, ignore_errors=True)
    if not moved:
        die(f"{key}: test mode produced no result CSV in {staging}")


# ---------------------------------------------------------------------------
# Training half (issue 04): sequential, one process at a time, resumable
# ---------------------------------------------------------------------------
def run_one(seed, mode, tag, device, pretrained, force):
    key = _split_key(seed, mode, tag)
    model_name = _model_name(key)
    weight = WEIGHTS_DIR / f"{model_name}_best.pth"
    result = RESULTS_DIR / f"{model_name}_train_result.csv"
    cfg = CONFIG_DIR / key  # without .yaml; load_config appends it
    staging = STAGING_ROOT / f"cm_{key}"

    if not force and weight.exists() and result.exists() and is_test_csv(result):
        log(f"Cached {key} — skipping")
        return
    if not force and weight.exists():
        log(f"{key}: weight cached, result missing — evaluating only")

    # --- training (skipped if the weight is cached) ---------------------
    if force or not weight.exists():
        shutil.rmtree(staging, ignore_errors=True)
        run_xpert("train", key, device, config=cfg,
                  extra=["--pretrained_model", pretrained,
                         "--num_epochs", str(NUM_EPOCHS)])
        if not staging.exists() or not list(staging.glob("*.pth")):
            die(f"{key}: training produced no weights in {staging}")
        WEIGHTS_DIR.mkdir(parents=True, exist_ok=True)
        items = (list(staging.glob(f"{model_name}_best.pth")) +
                 list(staging.glob(f"{model_name}_best.json")) +
                 list(staging.glob(f"{model_name}_fold_train_loss.png")))
        for item in items:
            dest = WEIGHTS_DIR / item.name
            if dest.exists():
                dest.unlink()
            shutil.move(str(item), str(dest))
        shutil.rmtree(staging, ignore_errors=True)
        log(f"{key}: trained -> {weight}")

    # --- test evaluation (skipped if a test-mode result CSV is cached) --
    if force or not (result.exists() and is_test_csv(result)):
        _eval_and_collect(
            key, staging, device, cfg,
            extra=["--saved_model_path", str(WEIGHTS_DIR),
                   "--pretrained_model", str(weight)])
        log(f"{key}: evaluated -> {result}")


def run_pretrain_anchor(seed, device, pretrained, force):
    key = f"ablation_pretrain_seed{seed}"
    model_name = _model_name(key)
    result = RESULTS_DIR / f"{model_name}_train_result.csv"
    cfg = CONFIG_DIR / key
    staging = STAGING_ROOT / f"cm_{key}"

    if not force and result.exists() and is_test_csv(result):
        log(f"Cached {key} (pretrain zero-shot) — skipping")
        return
    # --saved_model + --pretrained_model both = the step-04 L1000 pretrain
    # weight; the eval config mirrors the pretrain dims (240/21) so the
    # checkpoint loads tolerantly.
    _eval_and_collect(
        key, staging, device, cfg,
        extra=["--saved_model", pretrained,
               "--pretrained_model", pretrained])
    log(f"{key}: pretrain zero-shot evaluated -> {result}")


def run_training(masks_df, device, pretrained, force):
    WEIGHTS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    total = len(masks_df)
    for i, (_, row) in enumerate(masks_df.iterrows(), start=1):
        log(f"[{i}/{total}] {row['seed']} {row['mode']} "
            f"{row['fraction_tag']} (n_train={row['n_train']}, "
            f"n_test={row['n_test']})")
        run_one(int(row["seed"]), row["mode"], row["fraction_tag"], device,
                pretrained, force)
    for seed in SEEDS:
        run_pretrain_anchor(seed, device, pretrained, force)
    shutil.rmtree(WORKFLOW_ROOT / "finetune", ignore_errors=True)


# ---------------------------------------------------------------------------
# Analysis half (issue 05): aggregate -> plot -> report
# ---------------------------------------------------------------------------
def aggregate(masks_df):
    rows = []
    for _, mrow in masks_df.iterrows():
        key = _split_key(int(mrow["seed"]), mrow["mode"],
                         mrow["fraction_tag"])
        csv = RESULTS_DIR / f"{_model_name(key)}_train_result.csv"
        if not csv.exists():
            die(f"aggregate: missing result CSV {csv} — run training first")
        df = pd.read_csv(csv)
        if not is_test_csv(csv):
            die(f"aggregate: {csv} is not a test-mode result CSV — "
                f"re-run that (seed, mode, fraction)")
        met = df.iloc[0].to_dict()
        rows.append({
            "seed": int(mrow["seed"]),
            "mode": mrow["mode"],
            "fraction": float(mrow["fraction"]),
            "n_train": int(mrow["n_train"]),
            "n_test": int(mrow["n_test"]),
            "n_absent_drugs": int(mrow["n_absent_drugs"]),
            **{m: met.get(m, np.nan) for m in METRICS},
        })
    # pretrain anchors (fraction -> 0)
    for seed in SEEDS:
        key = f"ablation_pretrain_seed{seed}"
        csv = RESULTS_DIR / f"{_model_name(key)}_train_result.csv"
        if not csv.exists():
            die(f"aggregate: missing pretrain result CSV {csv}")
        df = pd.read_csv(csv)
        if not is_test_csv(csv):
            die(f"aggregate: {csv} is not a test-mode result CSV")
        met = df.iloc[0].to_dict()
        seed_test = masks_df[masks_df["seed"] == seed].iloc[0]
        rows.append({
            "seed": seed,
            "mode": "pretrain",
            "fraction": 0.0,
            "n_train": 0,
            "n_test": int(seed_test["n_test"]),
            "n_absent_drugs": int(seed_test["n_test_drugs"]),
            **{m: met.get(m, np.nan) for m in METRICS},
        })
    out = pd.DataFrame(rows, columns=[
        "seed", "mode", "fraction", "n_train", "n_test", "n_absent_drugs"] +
        METRICS)
    out = out.sort_values(["seed", "mode", "fraction"]).reset_index(drop=True)
    out.to_csv(OUTPUT_DIR / "ablation.csv", index=False)
    log(f"Aggregated {len(out)} rows -> {OUTPUT_DIR / 'ablation.csv'}")
    return out


def plot_ablation(ablation):
    anchor_x = min(FRACTIONS) / 2  # fraction -> 0 (below 1/16)
    x = np.array(FRACTIONS)
    fig, axes = plt.subplots(2, len(PLOT_METRICS), figsize=(19, 8))
    for r, mode in enumerate(MODES):
        sub = ablation[ablation["mode"] == mode]
        anchors = ablation[ablation["mode"] == "pretrain"]
        for c, (label, metric) in enumerate(PLOT_METRICS):
            ax = axes[r, c]
            means, stds = [], []
            for f in FRACTIONS:
                vals = sub.loc[sub["fraction"] == f, metric].astype(
                    float).values
                if len(vals) != len(SEEDS):
                    die(f"plot: expected {len(SEEDS)} seed points at "
                        f"{mode} {f} for {metric}, got {len(vals)}")
                means.append(vals.mean())
                stds.append(vals.std())
                ax.scatter([f] * len(vals), vals, s=18, color="black",
                           zorder=3, alpha=0.75)
            means, stds = np.array(means), np.array(stds)
            ax.plot(x, means, "-o", color="steelblue", lw=1.8, zorder=2)
            ax.fill_between(x, means - stds, means + stds, alpha=0.22,
                            color="steelblue")
            a_vals = anchors[metric].astype(float).values
            ax.scatter([anchor_x] * len(a_vals), a_vals, marker="D", s=42,
                       color="#e8a13c", zorder=4,
                       label=f"{ANCHOR_LABEL} (L1000, no finetune)")
            ax.axvline(anchor_x, color="#e8a13c", lw=0.8, ls=":", alpha=0.6)
            ax.set_xscale("log")
            ax.set_xticks(list(x) + [anchor_x])
            ax.set_xticklabels(
                ["1/16", "1/8", "1/4", "1/2", "1.0", ANCHOR_LABEL], fontsize=8)
            ax.set_title(f"{label} ({metric})", fontsize=11, fontweight="bold")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        axes[r, 0].set_ylabel(mode, fontsize=12, fontweight="bold")
        axes[r, 0].yaxis.set_label_position("left")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=1, fontsize=11,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("XPert finetune-corpus data-ablation — test metrics vs "
                 f"train-size fraction (3 seeds; x = fraction log scale; "
                 f"{ANCHOR_LABEL} = L1000 pretrain weight, no finetune, at "
                 f"fraction->0)", fontsize=13, fontweight="bold")
    fig.tight_layout()
    out_png = OUTPUT_DIR / "ablation.png"
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved: {out_png}")


def write_report(ablation):
    lines = [
        "# XPert finetune-corpus data-ablation — summary\n",
        "Diagnostic experiment (ADR 0018), decoupled from the official "
        "dtoxs-primary folds. 3 seeds x (warm-start | cold-drug) x "
        "5 train-size fractions {1/16, 1/8, 1/4, 1/2, 1.0}; each seed holds "
        "one unified test set (~15% of the 4-dataset union) shared by every "
        "fraction of both modes; train = a dataset-stratified f-subsample of "
        "the seed's train pool. Recipe identical to the official fold bench "
        "(4 epochs, lr 0.002, batch 32, cm.yaml, unimol+HG, step-04 L1000 "
        "pretrain init). The L1000-pretrain weight anchors the curve at 0 "
        "rows. This report changes no model used downstream.\n",
        f"Runs aggregated: {len(ablation)} (30 finetune + 3 pretrain "
        f"anchors).\n",
        "## Scale curve (mean of 3 seeds per mode)\n",
        "| mode | metric | f=1/16 | f=1/8 | f=1/4 | f=1/2 | f=1.0 | anchor |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for mode in MODES:
        sub = ablation[ablation["mode"] == mode]
        for label, metric in PLOT_METRICS:
            vals = [sub.loc[sub["fraction"] == f, metric].mean()
                    for f in FRACTIONS]
            anc = ablation[ablation["mode"] == "pretrain"][metric].mean()
            lines.append(f"| {mode} | {label} | " +
                         " | ".join(f"{v:.3f}" for v in vals) +
                         f" | {anc:.3f} |")
    lines += [
        "",
        "## Warm-start drug-coverage diagnostic (n_absent_drugs)",
        "",
        "Number of test drugs absent from the train subsample at that "
        "fraction (warm-start; cold-drug is strictly drug-disjoint by "
        "design, so its diagnostic always equals n_test_drugs).",
        "",
        "| seed | mode | fraction | n_train | n_test | n_absent_drugs |",
        "|---|---|---|---|---|---|",
    ]
    for _, r in ablation.iterrows():
        lines.append(f"| {r['seed']} | {r['mode']} | {r['fraction']} | "
                     f"{r['n_train']} | {r['n_test']} | "
                     f"{r['n_absent_drugs']} |")
    lines += [
        "",
        "## Interpretation guide",
        "",
        "- A metric that stops changing between f=1/2 and f=1.0 marks where "
        "more CM data stops paying off on the union-wide test distribution.",
        "- warm-start vs cold-drug spread shows whether scale helps known-drug "
        "(memorised) prediction differently from unseen-drug generalisation.",
        "- The anchor column is the un-finetuned L1000 pretrain weight "
        "(0 training rows), so the curve spans no-finetune to full-finetune.",
        "- A warm-start drop at small f may be 'forgot a known drug' — read "
        "the n_absent_drugs diagnostic before attributing it to scale.",
        "- **Warm vs cold scale ranges differ by design.** The cold-drug train "
        "pool strips every test drug's samples across all four datasets, so it "
        "holds only the non-test-drug rows (≈250 rows vs ≈3,700 warm): the "
        "cold curve's f=1.0 point is ~250 rows, NOT union-scale, and the two "
        "rows' x-ranges must not be read as comparable full-scale spans.",
    ]
    (OUTPUT_DIR / "ablation_report.md").write_text("\n".join(lines))
    log(f"Wrote {OUTPUT_DIR / 'ablation_report.md'}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Data-ablation benchmark for the v3 MIRI finetune "
                    "(ADR 0018)")
    parser.add_argument("--dry-run", action="store_true",
                        help="build masks + run the isolation gates, exit "
                             "without any GPU work")
    parser.add_argument("--pretrained-model", default=str(DEFAULT_PRETRAINED),
                        help="step-04 L1000 pretrain weight")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--force", action="store_true",
                        help="redo finished (seed, mode, fraction) runs / "
                             "regenerate the mask manifest")
    args = parser.parse_args()

    if not Path(str(FINETUNE_CONFIG)).exists():
        die(f"finetune config not found: {FINETUNE_CONFIG}")
    if not Path(str(PRETRAIN_EVAL_CONFIG)).exists():
        die(f"pretrain eval config not found: {PRETRAIN_EVAL_CONFIG}")
    if not UNION_H5AD.exists():
        die(f"union h5ad not found: {UNION_H5AD}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Mask construction reads only adata.obs — cheap, and byte-stable across
    # invocations with the same seeds.
    adata = sc.read_h5ad(UNION_H5AD, backed="r")
    log(f"Loaded {UNION_H5AD} (obs-only, {adata.n_obs} rows)")
    fresh = build_all_masks(adata)
    verify_gates(adata, fresh)
    masks_df = persist_masks(fresh, args.force)

    if args.dry_run:
        log("Dry-run complete — masks + gates validated, no GPU work "
            "performed.")
        return

    # Real run: full data for the temporary h5ads.
    adata = sc.read_h5ad(UNION_H5AD)
    prepare_staging(adata, masks_df)
    run_training(masks_df, args.device, args.pretrained_model, args.force)
    ablation = aggregate(masks_df)
    plot_ablation(ablation)
    write_report(ablation)
    log(f"Data-ablation complete — outputs in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
