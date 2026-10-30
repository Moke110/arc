#!/usr/bin/env python
"""Paired comparison of the 5k sensitivity run against the executed 18k run.

Reads only evaluated artifacts (`leaderboard.csv`, `per_drug.csv`); no model or
prediction is loaded.  The two runs share wells, folds and test drugs, so the
per-drug and per-fold differences are paired.

Usage
-----
    python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/compare_18k.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from repo_root import BENCH  # noqa: E402

DEFAULT_18K = BENCH.parent / "03.1_cold_drug_18k" / "output" / "results"


def _leaderboard(path: Path, prefix: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    keep = ["model", "headline_auprc", "mean_k100", "mean_k500", "pcc", "spearman", "r2", "rmse", "rank"]
    return frame[keep].rename(columns={c: f"{prefix}_{c}" for c in keep if c != "model"})


def compare(fivek: Path, eighteenk: Path) -> pd.DataFrame:
    five = _leaderboard(fivek / "leaderboard.csv", "k5")
    eight = _leaderboard(eighteenk / "leaderboard.csv", "k18")
    merged = five.merge(eight, on="model", how="outer", validate="one_to_one")
    merged["delta_headline"] = merged["k5_headline_auprc"] - merged["k18_headline_auprc"]
    merged["delta_rank"] = merged["k5_rank"] - merged["k18_rank"]

    def per_drug(root: Path, prefix: str) -> pd.DataFrame:
        frame = pd.read_csv(root / "per_drug.csv")
        grouped = frame.groupby(["model", "drug_key"], observed=True)["auprc"].mean().reset_index()
        return grouped.rename(columns={"auprc": f"{prefix}_auprc"})

    drugs = per_drug(fivek, "k5").merge(
        per_drug(eighteenk, "k18"), on=["model", "drug_key"], how="inner", validate="one_to_one"
    )
    drugs["delta"] = drugs["k5_auprc"] - drugs["k18_auprc"]

    rng = np.random.default_rng(20261022)
    rows = []
    for model, group in drugs.groupby("model", observed=True):
        values = group["delta"].to_numpy(dtype=np.float64)
        draws = rng.choice(values, size=(2000, len(values)), replace=True).mean(axis=1)
        rows.append(
            {
                "model": model,
                "mean_drug_delta": float(values.mean()),
                "median_drug_delta": float(np.median(values)),
                "n_drugs": int(len(values)),
                "share_drugs_up": float((values > 0).mean()),
                "delta_ci_low": float(np.quantile(draws, 0.025)),
                "delta_ci_high": float(np.quantile(draws, 0.975)),
            }
        )
    summary = pd.DataFrame(rows)
    return merged.merge(summary, on="model", how="left").sort_values("k5_rank")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fivek", type=Path, default=BENCH / "output" / "results")
    parser.add_argument("--eighteenk", type=Path, default=DEFAULT_18K)
    parser.add_argument("--out", type=Path, default=BENCH / "output" / "results" / "comparison_18k.csv")
    args = parser.parse_args(argv)
    table = compare(args.fivek, args.eighteenk)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, index=False)
    columns = [
        "model", "k5_rank", "k18_rank", "k5_headline_auprc", "k18_headline_auprc",
        "delta_headline", "mean_drug_delta", "delta_ci_low", "delta_ci_high", "share_drugs_up",
    ]
    print(table[columns].to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
