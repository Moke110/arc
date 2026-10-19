#!/usr/bin/env python3
"""Step 05 — Data-ablation benchmark for the v2 MIRI finetune (ADR 0019).

Measures how the finetune-corpus SCALE affects XPert test performance in
warm-start and cold-drug modes, **decoupled from the official folds** (the
`full/` bench, ADR 0018 term *dtoxs-primary benchmark*). The experiment is a
diagnostic analysis only: it changes no model used downstream (step 06 cm_full
and steps 07-11 are untouched).

v2's union is small — 958 rows / 45 unique drugs across two datasets
(dtoxs_primary 662 + dtoxs_ipsc 296) — so a row-uniform test (v3's recipe)
would draw ~40 of 45 drugs into test and collapse the cold-drug train pool to
a handful of rows. The v2 design therefore **partitions by drug** instead:

Design (3 seeds x 2 modes x 5 fractions = 30 finetune runs + 3 anchors):
  - The 45 union drugs are split into a **disjoint tripartition** (3 x 15
    test drugs, fixed master RNG): seed s holds drug-third s as its *test
    drugs*, so every union drug is tested exactly once across the 3 seeds and
    the seeds' test sets never overlap.
  - Per seed, the **unified test set** (~15% of the union, drawn only from
    that seed's 15 test drugs, per-drug allocation proportional to each test
    drug's row count minus one) is shared by every fraction of both modes.
    Each test drug keeps >= 1 row in the pool (warm-start f=1 overlap 1.0).
  - Train pools: warm-start = union minus test (all 45 drugs, ~814 rows);
    cold-drug = union minus test minus every row of the 15 test drugs (the 30
    hold-out drugs only, ~610-655 rows, strictly drug-disjoint at every
    fraction).
  - Fractions {20%, 40%, 60%, 80%, 100%}: each a dataset-stratified
    f-subsample of the seed's pool (per-dataset row counts x f, composition
    held constant), drawn with a derived deterministic seed per
    (seed, mode, fraction). 100% is re-run on this scheme (never the full/
    bench).
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
  profiles/            per-run per-cell prediction profiles (*_predict_profile.npy)
  ablation.csv         33 rows (30 runs + 3 pretrain anchors) + VCC columns
  ablation.png         4 rows x 4 cols (rows 0/2 = warm-start|cold-drug deg
                       metrics RMSE|PCC|Spearman|R2; rows 1/3 = the same
                       modes' VCC metrics PDS|JAC|FID|REACH; x = fraction on a
                       log scale; 3 seed points + mean +/- std; pretrain
                       anchor at fraction -> 0)
  ablation_report.md   scale-curve + per-fraction drug-coverage summary

Per-run finetune YAMLs and the pretrain-eval config both live at
05_xpert_finetune_bench/config/ (exactly v3's nesting: the XPert loader
resolves a config's relative dataset paths against its third parent, config
-> stage -> workflow root; the shared cm_pretrain_eval.yaml is consumed by
both the full/ bench and the ablation anchors, no copy).

Usage:
    python 05_xpert_finetune_bench/data-ablation/xpert_data_ablation.py --dry-run
    python 05_xpert_finetune_bench/data-ablation/xpert_data_ablation.py
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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "_lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "full"))
from utils import WORKFLOW_ROOT  # noqa: E402
from xpert_finetune_bench import (  # noqa: E402
    is_test_csv,
    run_xpert,
    STAGING_ROOT,
    DEFAULT_PRETRAINED,
)
from vcc_metrics import compute_vcc_metrics  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "05_xpert_finetune_bench"
OUTPUT_DIR = STEP_DIR / "data-ablation"
WORK_DIR = OUTPUT_DIR / "work"
WEIGHTS_DIR = OUTPUT_DIR / "weights"
RESULTS_DIR = OUTPUT_DIR / "results"
PROFILES_DIR = OUTPUT_DIR / "profiles"
# Generated per-run YAMLs live at 05_xpert_finetune_bench/config/ (v3's
# nesting): the XPert loader resolves a config's relative dataset paths
# (DTOX_data_root, HG/UniMol/SMILES drug embeds, PPI vector) against its third
# parent (config -> stage -> workflow root), so the YAMLs must sit two levels
# under the workflow root.
CONFIG_DIR = STEP_DIR / "config"
UNION_H5AD = WORKFLOW_ROOT / "03_cm_finetune_data" / "output" / "cm_union.h5ad"
FINETUNE_CONFIG = WORKFLOW_ROOT / "06_xpert_finetune" / "config" / "cm.yaml"
# The pretrain-eval config is shared with the full/ fold bench (no copy).
# It stays at the step-top-level config/ for the loader's third-parent rule.
PRETRAIN_EVAL_CONFIG = STEP_DIR / "config" / "cm_pretrain_eval.yaml"

SEEDS = [1, 2, 3]
MODES = ["warm-start", "cold-drug"]
FRACTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
FRAC_TAGS = ["20", "40", "60", "80", "100"]
TEST_FRACTION = 0.15
TEST_FRACTION_TOL = 0.01
NUM_EPOCHS = 4
# Fixed master RNG for the drug tripartition: independent of the per-run
# (seed, mode, fraction) draws, so the test-drug thirds are stable across
# invocations and only the subsample draws consume the seed-derived RNGs.
TRIPARTITION_RNG_SEED = 0
MASK_CSV = OUTPUT_DIR / "ablation_masks.csv"
METRICS = ["MSE_deg", "RMSE_deg", "MAE_deg", "Pearson_deg", "Spearman_deg",
           "R2_deg"]
# (panel label, metric column) — the four headline panels per mode.
PLOT_METRICS = [("RMSE", "RMSE_deg"), ("PCC", "Pearson_deg"),
                ("Spearman", "Spearman_deg"), ("R2", "R2_deg")]
# VCC 2026 metrics (raw, higher-is-better): PDS / JAC / FID / REACH.
VCC_METRICS = ["PDS", "JAC", "FID", "REACH"]
PLOT_VCC = [(m, m) for m in VCC_METRICS]
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
def build_tripartition(pert):
    """Disjoint tripartition of the union's drugs into SEEDS equal thirds.

    Deterministic under TRIPARTITION_RNG_SEED: shuffle the sorted drug list
    once, then chunk into len(SEEDS) equal groups; seed i takes group i. Every
    union drug lands in exactly one seed's test-drug set.
    """
    drugs = np.array(sorted(set(int(p) for p in pert)))
    if len(drugs) % len(SEEDS) != 0:
        die(f"cannot tripartition {len(drugs)} drugs into {len(SEEDS)} equal "
            f"thirds")
    rng = np.random.default_rng(TRIPARTITION_RNG_SEED)
    shuffled = rng.permutation(drugs)
    groups = np.array_split(shuffled, len(SEEDS))
    sizes = {len(g) for g in groups}
    if len(sizes) != 1:
        die(f"uneven tripartition: group sizes {sorted(len(g) for g in groups)}")
    return {seed: set(int(d) for d in groups[i]) for i, seed in enumerate(SEEDS)}


def derived_seed(seed, mode, fraction):
    """Deterministic derived RNG seed per (seed, mode, fraction)."""
    h = hashlib.md5(f"{seed}|{mode}|{fraction}".encode()).hexdigest()
    return int(h, 16)


def build_test_indices(adata, test_drugs, rng):
    """Unified test set (~15% of the union) drawn from the seed's test drugs.

    Per-drug allocation: each test drug contributes >= 1 test cell, and the
    remaining budget is apportioned over the test drugs proportional to each
    drug's *eligible* rows (row count minus one — one row is always left in
    the pool so the warm-start f=1 train pool keeps every test drug, overlap
    1.0). Largest remainder (Hamilton) over the excess capacity, total exactly
    round(TEST_FRACTION * n) (or the eligible total, whichever is smaller).
    The specific rows within each drug are chosen by ``rng`` (deterministic
    per seed).
    """
    pert = adata.obs["pert_idx"].astype(int).values
    n = adata.n_obs
    by_drug = defaultdict(list)
    for i in range(n):
        by_drug[int(pert[i])].append(i)
    drugs = sorted(test_drugs)
    eligible = {d: len(by_drug[d]) - 1 for d in drugs}
    total_budget = int(round(TEST_FRACTION * n))
    base = 1  # one guaranteed test cell per test drug
    n_drugs = len(drugs)
    if sum(eligible.values()) < n_drugs:
        die(f"test drugs have no spare rows: {sum(eligible.values())} "
            f"eligible < {n_drugs} drugs")
    extra_budget = max(0, min(total_budget - n_drugs,
                              sum(eligible.values()) - n_drugs))
    # weights for the excess seats = eligible minus the base seat
    weights = {d: eligible[d] - base for d in drugs}
    extra = {d: 0 for d in drugs}
    if extra_budget > 0:
        wsum = sum(weights.values())
        exact = {d: weights[d] * extra_budget / wsum for d in drugs}
        extra = {d: int(exact[d]) for d in drugs}
        rem = {d: exact[d] - extra[d] for d in drugs}
        shortfall = extra_budget - sum(extra.values())
        for d in sorted(drugs, key=lambda d: (-rem[d], d)):
            if shortfall <= 0:
                break
            if extra[d] < weights[d]:
                extra[d] += 1
                shortfall -= 1
    test = set()
    for d in drugs:
        k = base + extra[d]
        if k <= 0:
            continue
        pool = np.array(sorted(by_drug[d]), dtype=int)
        chosen = rng.choice(pool, size=min(k, len(pool)), replace=False)
        test.update(chosen.tolist())
    return frozenset(test)


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


def build_all_masks(adata, tripartition):
    """Per-(seed, mode, fraction) masks for all seeds. Deterministic."""
    pert = adata.obs["pert_idx"].astype(int).values
    ds = adata.obs["dataset_id"].astype(str).values
    wells = adata.obs["well"].astype(str).values
    n = adata.n_obs
    idx = np.arange(n)

    rows = []
    for seed in SEEDS:
        test_drugs = tripartition[seed]
        rng = np.random.default_rng(seed)
        test = build_test_indices(adata, test_drugs, rng)
        in_test = np.zeros(n, dtype=bool)
        in_test[list(test)] = True
        warm_pool = idx[~in_test]
        keep = ~np.isin(pert, np.array(sorted(test_drugs)))
        cold_pool = idx[keep]
        for mode in MODES:
            pool = warm_pool if mode == "warm-start" else cold_pool
            for frac, tag in zip(FRACTIONS, FRAC_TAGS):
                rng2 = np.random.default_rng(derived_seed(seed, mode, frac))
                train = _stratify_sub(pool, ds, frac, rng2)
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


def verify_gates(adata, masks_df, tripartition):
    """Assert the ADR 0019 design invariants; exit 1 on any failure."""
    pert = adata.obs["pert_idx"].astype(int).values
    ds = adata.obs["dataset_id"].astype(str).values
    well_to_idx = {w: i for i, w in
                   enumerate(adata.obs["well"].astype(str).values)}
    n = adata.n_obs
    problems = []

    # tripartition: every union drug in exactly one seed's test-drug set
    all_drugs = set(int(p) for p in pert)
    per_seed = {seed: tripartition[seed] for seed in SEEDS}
    covered = set().union(*per_seed.values())
    if covered != all_drugs:
        problems.append(f"tripartition coverage {len(covered)} != "
                        f"{len(all_drugs)} union drugs")
    if len(all_drugs) != sum(len(s) for s in per_seed.values()):
        problems.append("tripartition sets overlap across seeds")

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
        if test_drugs != tripartition[seed]:
            problems.append(f"seed {seed}: test wells carry drugs "
                            f"{sorted(test_drugs - tripartition[seed])[:5]} "
                            f"outside the seed's test-drug third")

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

        # every test drug keeps >=1 test cell and >=1 warm-pool row
        per_drug_test = defaultdict(int)
        for w in test_wells:
            per_drug_test[int(pert[well_to_idx[w]])] += 1
        for d in tripartition[seed]:
            if per_drug_test.get(d, 0) < 1:
                problems.append(f"seed {seed}: test drug {d} has 0 test cells")
            warm_rows_of_d = sum(
                1 for i in np.where(pert == d)[0] if not in_test_mask[i])
            if warm_rows_of_d < 1:
                problems.append(f"seed {seed}: test drug {d} has 0 warm-pool "
                                f"rows (f=1 overlap broken)")

    if problems:
        log("SPLIT GATE FAILED:", "ERROR")
        for p in problems:
            log(f"  - {p}", "ERROR")
        sys.exit(1)
    log(f"Split gates passed for {len(SEEDS)} seeds: drug tripartition "
        f"covers all {len(all_drugs)} union drugs exactly once; test size "
        f"≈{TEST_FRACTION}; train/test well-disjoint; warm-start f=1 overlap "
        f"1.0; cold-drug overlap 0.0 at every fraction; dataset-stratified "
        f"composition stable within the apportionment bound (1/n_train).")
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
    # string columns round-trip through CSV type inference (fraction_tag may
    # come back numeric) — compare on canonical string form.
    for c in ("fraction_tag", "test_wells", "train_wells"):
        a[c] = a[c].astype(str)
        b[c] = b[c].astype(str)
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
    RESULTS_DIR and the per-cell prediction profile to PROFILES_DIR. Shared by
    the finetune runs and the pretrain anchors."""
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True, exist_ok=True)
    run_xpert("test", key, device, config=config, extra=extra)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PROFILES_DIR.mkdir(parents=True, exist_ok=True)
    produced = list(staging.glob("*_train_result.csv"))
    if not produced:
        die(f"{key}: test mode produced no result CSV in {staging}")
    for csv in produced:
        dest = RESULTS_DIR / csv.name
        if dest.exists():
            csv.unlink()
        else:
            shutil.move(str(csv), str(dest))
    for npy in staging.glob("*_predict_profile.npy"):
        dest = PROFILES_DIR / f"{key}_predict_profile.npy"
        if dest.exists():
            dest.unlink()
        shutil.move(str(npy), str(dest))
    shutil.rmtree(staging, ignore_errors=True)


# ---------------------------------------------------------------------------
# Training half: sequential, one process at a time, resumable
# ---------------------------------------------------------------------------
def run_one(seed, mode, tag, device, pretrained, force):
    key = _split_key(seed, mode, tag)
    model_name = _model_name(key)
    weight = WEIGHTS_DIR / f"{model_name}_best.pth"
    result = RESULTS_DIR / f"{model_name}_train_result.csv"
    profile = PROFILES_DIR / f"{key}_predict_profile.npy"
    cfg = CONFIG_DIR / key  # without .yaml; load_config appends it
    staging = STAGING_ROOT / f"cm_{key}"

    if (not force and weight.exists() and result.exists()
            and is_test_csv(result) and profile.exists()):
        log(f"Cached {key} — skipping")
        return
    if not force and weight.exists():
        log(f"{key}: weight cached, result/profile missing — evaluating only")

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

    # --- test evaluation (skipped if result CSV + prediction profile are
    # both cached) --------------------------------------------------------
    if force or not (result.exists() and is_test_csv(result)) or \
            not profile.exists():
        _eval_and_collect(
            key, staging, device, cfg,
            extra=["--saved_model_path", str(WEIGHTS_DIR),
                   "--pretrained_model", str(weight),
                   "--output_profile", "True"])
        log(f"{key}: evaluated -> {result}")


def run_pretrain_anchor(seed, device, pretrained, force):
    key = f"ablation_pretrain_seed{seed}"
    model_name = _model_name(key)
    result = RESULTS_DIR / f"{model_name}_train_result.csv"
    profile = PROFILES_DIR / f"{key}_predict_profile.npy"
    cfg = CONFIG_DIR / key
    staging = STAGING_ROOT / f"cm_{key}"

    if (not force and result.exists() and is_test_csv(result)
            and profile.exists()):
        log(f"Cached {key} (pretrain zero-shot) — skipping")
        return
    # --saved_model + --pretrained_model both = the step-04 L1000 pretrain
    # weight; the eval config mirrors the pretrain dims (240/21) so the
    # checkpoint loads tolerantly.
    _eval_and_collect(
        key, staging, device, cfg,
        extra=["--saved_model", pretrained,
               "--pretrained_model", pretrained,
               "--output_profile", "True"])
    log(f"{key}: pretrain zero-shot evaluated -> {result}")


def run_training(masks_df, device, pretrained, force):
    WEIGHTS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PROFILES_DIR.mkdir(parents=True, exist_ok=True)
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
# Analysis half: aggregate -> plot -> report
# ---------------------------------------------------------------------------
def _vcc_for_run(seed, mode, tag):
    """Compute the four raw VCC metrics for one (seed, mode, fraction) run or
    pretrain anchor from its cached per-cell prediction profile, grouped by
    drug (pert_idx)."""
    if mode == "pretrain":
        key = f"ablation_pretrain_seed{seed}"
        temp = WORK_DIR / f"seed{seed}_pretrain.h5ad"
        fold = "ablation_test"
    else:
        key = _split_key(seed, mode, tag)
        temp = WORK_DIR / f"seed{seed}_{mode}.h5ad"
        fold = f"ablation_{tag}"
    profile = PROFILES_DIR / f"{key}_predict_profile.npy"
    if not profile.exists():
        die(f"aggregate: missing prediction profile {profile} — run "
            f"inference first")
    prof = np.load(profile, allow_pickle=True).item()
    deg_true = np.asarray(prof["deg_true"], dtype=np.float64)
    deg_pred = np.asarray(prof["deg_pred"], dtype=np.float64)
    adata = sc.read_h5ad(temp, backed="r")
    drugs = adata.obs.loc[adata.obs[fold] == "test", "pert_idx"].astype(
        int).values
    if len(drugs) != deg_true.shape[0]:
        die(f"aggregate: profile/cell mismatch for {key}: "
            f"{deg_true.shape[0]} predictions vs {len(drugs)} test cells")
    return compute_vcc_metrics(deg_true, deg_pred, drugs)


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
        vcc = _vcc_for_run(int(mrow["seed"]), mrow["mode"],
                           mrow["fraction_tag"])
        rows.append({
            "seed": int(mrow["seed"]),
            "mode": mrow["mode"],
            "fraction": float(mrow["fraction"]),
            "n_train": int(mrow["n_train"]),
            "n_test": int(mrow["n_test"]),
            "n_absent_drugs": int(mrow["n_absent_drugs"]),
            **{m: met.get(m, np.nan) for m in METRICS},
            **{m: vcc.get(m, np.nan) for m in VCC_METRICS},
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
        vcc = _vcc_for_run(seed, "pretrain", None)
        rows.append({
            "seed": seed,
            "mode": "pretrain",
            "fraction": 0.0,
            "n_train": 0,
            "n_test": int(seed_test["n_test"]),
            "n_absent_drugs": int(seed_test["n_test_drugs"]),
            **{m: met.get(m, np.nan) for m in METRICS},
            **{m: vcc.get(m, np.nan) for m in VCC_METRICS},
        })
    out = pd.DataFrame(rows, columns=[
        "seed", "mode", "fraction", "n_train", "n_test", "n_absent_drugs"] +
        METRICS + VCC_METRICS)
    out = out.sort_values(["seed", "mode", "fraction"]).reset_index(drop=True)
    out.to_csv(OUTPUT_DIR / "ablation.csv", index=False)
    log(f"Aggregated {len(out)} rows -> {OUTPUT_DIR / 'ablation.csv'}")
    return out


def plot_ablation(ablation):
    # The L1000-pretrain anchor is plotted as a regular sixth column at
    # fraction -> 0 (x = min(FRACTIONS)/2, since 0 sits off a log axis),
    # tick-labelled "0" and styled exactly like the other columns: black
    # per-seed scatter + blue mean +/- std curve, no separate legend entry.
    # Four rows: warm-start (rows 1-2) then cold-drug (rows 3-4); within each
    # mode, the original deg metrics row is followed by the VCC metrics row.
    anchor_x = min(FRACTIONS) / 2
    x = np.array(FRACTIONS)
    x_plot = np.concatenate(([anchor_x], x))
    x_labels = ["0", "20%", "40%", "60%", "80%", "100%"]
    rows = [(m, metric_list)
            for m in MODES
            for metric_list in (PLOT_METRICS, PLOT_VCC)]
    fig, axes = plt.subplots(len(rows), len(PLOT_METRICS), figsize=(19, 16))
    for r, (mode, metric_list) in enumerate(rows):
        sub = ablation[ablation["mode"] == mode]
        anchors = ablation[ablation["mode"] == "pretrain"]
        for c, (label, metric) in enumerate(metric_list):
            ax = axes[r, c]
            means, stds = [], []
            for j, f in enumerate(x_plot):
                if j == 0:
                    vals = anchors[metric].astype(float).values
                else:
                    vals = sub.loc[sub["fraction"] == x[j - 1], metric].astype(
                        float).values
                if len(vals) != len(SEEDS):
                    die(f"plot: expected {len(SEEDS)} seed points at "
                        f"{mode} {f} for {metric}, got {len(vals)}")
                means.append(vals.mean())
                stds.append(vals.std())
                ax.scatter([f] * len(vals), vals, s=18, color="black",
                           zorder=3, alpha=0.75)
            means, stds = np.array(means), np.array(stds)
            ax.plot(x_plot, means, "-o", color="steelblue", lw=1.8, zorder=2)
            ax.fill_between(x_plot, means - stds, means + stds, alpha=0.22,
                            color="steelblue")
            ax.set_xscale("log")
            ax.set_xticks(list(x_plot))
            ax.set_xticklabels(x_labels, fontsize=8)
            ax.set_title(f"{label} ({metric})", fontsize=11, fontweight="bold")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        axes[r, 0].set_ylabel(mode, fontsize=12, fontweight="bold")
        axes[r, 0].yaxis.set_label_position("left")
    fig.suptitle("XPert finetune-corpus data-ablation — test metrics vs "
                 "train-size fraction (3 seeds; x = fraction log scale; "
                 "column 0 = L1000 pretrain weight, no finetune; VCC rows = "
                 "PDS/JAC/FID/REACH, raw, higher-is-better)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    out_png = OUTPUT_DIR / "ablation.png"
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved: {out_png}")


def write_report(ablation):
    lines = [
        "# XPert finetune-corpus data-ablation — summary\n",
        "Diagnostic experiment (ADR 0019), decoupled from the official "
        "full/ folds. v2's union is 958 rows / 45 drugs across the two "
        "finetune datasets (dtoxs_primary + dtoxs_ipsc), so the per-seed "
        "test sets are **drug-partitioned** (not row-uniform as in v3's "
        "4-dataset ablation): the 45 drugs are split into a disjoint "
        "tripartition of 3 x 15 test drugs; each seed tests one third "
        "(every union drug tested exactly once across the 3 seeds), with a "
        "~15%-of-union unified test set drawn only from that third. Train = a "
        "dataset-stratified f-subsample of the seed's pool over 5 train-size "
        "fractions {20%, 40%, 60%, 80%, 100%}; warm-start = union minus test "
        "(~814 rows, all 45 drugs); cold-drug = the 30 hold-out drugs only "
        "(~610-655 rows, strictly drug-disjoint at every fraction). Recipe "
        "identical to the official fold bench (4 epochs, lr 0.002, batch 32, "
        "cm.yaml, unimol+HG, step-04 L1000 pretrain init). The L1000-pretrain "
        "weight anchors the curve at 0 rows. This report changes no model "
        "used downstream.\n",
        f"Runs aggregated: {len(ablation)} (30 finetune + 3 pretrain "
        f"anchors).\n",
        "## Scale curve (mean of 3 seeds per mode)\n",
        "| mode | metric | f=20% | f=40% | f=60% | f=80% | f=100% | anchor |",
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
        for label, metric in PLOT_VCC:
            vals = [sub.loc[sub["fraction"] == f, metric].mean()
                    for f in FRACTIONS]
            anc = ablation[ablation["mode"] == "pretrain"][metric].mean()
            lines.append(f"| {mode} | {label} | " +
                         " | ".join(f"{v:.3f}" for v in vals) +
                         f" | {anc:.3f} |")
    lines += [
        "",
        "PDS / JAC / FID / REACH are the four raw (unscaled) higher-is-better "
        "Virtual Cell Challenge 2026 metrics (cell_eval2 vcc2026 profile), "
        "computed per drug from the per-cell predicted DEGs and averaged over "
        "the test drugs. Because the test set holds ~10 cells per test drug, "
        "the DE gene set for JAC/FID/REACH uses the top-50 genes by |log2FC| "
        "(not a Wilcoxon p_adj<0.05 set, which is underpowered at n=10); PDS "
        "is a cosine-retrieval rank and needs no significance test. No "
        "target-gene exclusion applies (drug perturbations). Drugs with very "
        "few rows (< 10) contribute proportionally fewer test cells.",
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
        "- A metric that stops changing between f=80% and f=100% marks where "
        "more CM data stops paying off on the union-wide test distribution.",
        "- warm-start vs cold-drug spread shows whether scale helps known-drug "
        "(memorised) prediction differently from unseen-drug generalisation.",
        "- The anchor column is the un-finetuned L1000 pretrain weight "
        "(0 training rows), so the curve spans no-finetune to full-finetune.",
        "- A warm-start drop at small f may be 'forgot a known drug' — read "
        "the n_absent_drugs diagnostic before attributing it to scale.",
        "- **Warm vs cold scale ranges differ by design.** The cold-drug "
        "train pool strips every row of the 15 test drugs, so it holds only "
        "the 30 hold-out drugs (~610-655 rows vs ~814 warm): the cold curve's "
        "f=100% point is the hold-out-drug scale, NOT union scale, and the "
        "two rows' x-ranges must not be read as comparable full-scale spans.",
    ]
    (OUTPUT_DIR / "ablation_report.md").write_text("\n".join(lines))
    log(f"Wrote {OUTPUT_DIR / 'ablation_report.md'}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Data-ablation benchmark for the v2 MIRI finetune "
                    "(ADR 0019)")
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
    pert = adata.obs["pert_idx"].astype(int).values
    tripartition = build_tripartition(pert)
    for seed in SEEDS:
        log(f"seed {seed}: {len(tripartition[seed])} test drugs "
            f"(of {len(set(int(p) for p in pert))} union drugs)")
    fresh = build_all_masks(adata, tripartition)
    verify_gates(adata, fresh, tripartition)
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
