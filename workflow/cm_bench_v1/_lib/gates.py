"""
Gate assertions for the Stage 2 v1 benchmark.

The repository has no unit-test framework; correctness is enforced by gates that
fail loudly (ticket conventions). These are the invariant checks the benchmark
data build, splits and metric must satisfy before any result is trusted.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from metric import direction_aware_auprc, prevalence  # noqa: E402


def gate(cond, msg):
    if not cond:
        raise AssertionError(msg)


def check_baselines_finite(baseline):
    gate(np.isfinite(baseline).all(),
         "per-plate baseline has non-finite entries (treated well without a "
         "usable baseline)")


def check_universe(universe):
    gate(len(universe) > 0, "empty common gene universe")
    gate(len(universe) == len(set(universe)), "duplicate gene symbols in universe")


def check_eligibility(groups, fold_df, min_wells):
    """Every drug is recorded; eligibility matches the >= min-wells rule."""
    per = groups.groupby(["dataset", "drug"])["n_wells"].sum()
    elig = set(per[per >= min_wells].index.get_level_values("drug"))
    gate(len(elig) > 0, "no eligible drugs (>= %d wells) found" % min_wells)
    gate(set(fold_df.loc[fold_df["eligible"], "drug"]) == elig,
         "eligibility in the split mask disagrees with the >= min-wells rule")
    gate(set(fold_df["drug"]) == set(groups["drug"]),
         "fold table does not record every union drug")


def check_cold_drug_masks(fold_df, groups):
    elig = fold_df[fold_df["eligible"]]
    gate((elig["fold"] >= 0).all(), "eligible drug without a fold")
    gate(elig["fold"].nunique() > 1, "cold-drug split has a single fold")
    gate(not elig["drug"].duplicated().any(), "drug assigned to >1 fold")
    # zero train/test drug overlap in every fold
    folds = sorted(elig["fold"].unique())
    for f in folds:
        held = set(fold_df.loc[(fold_df["fold"] == f) & fold_df["eligible"], "drug"])
        test = set(groups.loc[groups["drug"].isin(held), "drug"])
        train = set(groups.loc[~groups["drug"].isin(held), "drug"])
        gate(not (test & train), f"fold {f}: train/test drug overlap")
    # every eligible drug tested exactly once
    covered = []
    for f in folds:
        covered += list(fold_df.loc[(fold_df["fold"] == f) & fold_df["eligible"], "drug"])
    gate(sorted(covered) == sorted(elig["drug"]),
         "cold-drug folds do not cover each eligible drug exactly once")


def check_cold_dose_masks(dose_df, groups):
    if dose_df.empty:
        return
    for split in ("interpolation", "extrapolation"):
        held = set(dose_df.loc[dose_df["held_out"] & (dose_df["split"] == split),
                               "dose_level"])
        gate(len(held) == 1, f"cold-dose {split}: expected exactly one held-out level")
        train_levels = set(dose_df.loc[~dose_df["held_out"]
                                       & (dose_df["split"] == split), "dose_level"])
        gate(not (held & train_levels),
             f"cold-dose {split}: held-out level present in train")


def check_metric():
    rng = np.random.default_rng(0)
    n = 2000
    labels = np.zeros(n, dtype=np.int8)
    up_idx = rng.choice(n, 60, replace=False)
    labels[up_idx] = 1
    down_idx = rng.choice(np.setdiff1d(np.arange(n), up_idx), 60, replace=False)
    labels[down_idx] = -1
    # perfect ranking: up genes highest, down genes lowest
    perfect = labels.astype(float) + 1e-3 * rng.standard_normal(n)
    gate(abs(direction_aware_auprc(perfect, labels) - 1.0) < 1e-9,
         "perfect ranking must give AUPRC 1.0")
    zeros = np.zeros(n)
    prev = prevalence(labels)
    gate(abs(direction_aware_auprc(zeros, labels) - prev) < 1e-9,
         "all-zeros prediction must give the prevalence baseline")
    rand = rng.standard_normal(n)
    a = direction_aware_auprc(rand, labels)
    gate(abs(a - prev) < 0.05,
         f"random ranking AUPRC {a:.3f} far from prevalence {prev:.3f}")
