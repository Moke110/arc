"""
Persist / load the Stage 2 v1 benchmark artifacts.

Layout (all under `output/bench/`):

    samples.csv          one row per treated well (dataset, drug, dose, ...)
    universe.txt         common gene universe, one symbol per line
    matrices.npz         baseline (n, G) float32, delta (n, G) float32
    groups.csv           one row per (dataset, drug, dose) group
    true_sets.npz        labels (n_groups, G) int8  (1 up, -1 down, 0 n.s.)
    splits/cold_drug.csv five cold-drug folds
    splits/cold_dose.csv leave-dose-level-out masks
    dataset_summary.csv  per-dataset ingestion stats
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import BENCH_DIR  # noqa: E402

SAMPLE_COLS = ["well", "dataset", "drug", "compound_code", "smiles", "dose_uM",
               "log10_dose", "dose_level", "plate", "time", "group_key"]


def _ensure():
    BENCH_DIR.mkdir(parents=True, exist_ok=True)
    (BENCH_DIR / "splits").mkdir(parents=True, exist_ok=True)


def save_bench(meta, universe, baseline, delta):
    _ensure()
    meta.to_csv(BENCH_DIR / "samples.csv", index=False)
    (BENCH_DIR / "universe.txt").write_text("\n".join(universe) + "\n")
    np.savez_compressed(BENCH_DIR / "matrices.npz",
                        baseline=baseline.astype(np.float32),
                        delta=delta.astype(np.float32))


def load_bench():
    meta = pd.read_csv(BENCH_DIR / "samples.csv", keep_default_na=False, na_values=[""])
    universe = (BENCH_DIR / "universe.txt").read_text().split()
    z = np.load(BENCH_DIR / "matrices.npz")
    return {"meta": meta, "universe": universe,
            "baseline": z["baseline"].astype(float),
            "delta": z["delta"].astype(float)}


def save_true_sets(groups, labels):
    _ensure()
    groups.to_csv(BENCH_DIR / "groups.csv", index=False)
    np.savez_compressed(BENCH_DIR / "true_sets.npz", labels=labels.astype(np.int8))


def load_true_sets():
    groups = pd.read_csv(BENCH_DIR / "groups.csv")
    labels = np.load(BENCH_DIR / "true_sets.npz")["labels"]
    return groups, labels


def save_cold_drug(fold_df):
    _ensure()
    fold_df.to_csv(BENCH_DIR / "splits" / "cold_drug.csv", index=False)


def save_cold_dose(dose_df):
    _ensure()
    dose_df.to_csv(BENCH_DIR / "splits" / "cold_dose.csv", index=False)


def save_dataset_summary(df):
    _ensure()
    df.to_csv(BENCH_DIR / "dataset_summary.csv", index=False)


def group_key(dataset, drug, dose_level):
    return f"{dataset}|{drug}|{int(dose_level)}"
