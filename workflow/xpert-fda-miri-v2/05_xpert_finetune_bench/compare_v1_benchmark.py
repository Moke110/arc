#!/usr/bin/env python3
"""
Benchmark comparison: v2 union-CM warm-start vs v1 DToxS warm-start.

Writes 05_xpert_finetune_bench/output/benchmark_vs_v1.csv — the mean ± std per
deg metric for the union-CM (v2) and DToxS-only (v1) warm-start benchmarks,
plus the per-fold tables, so the union finetune's effect on prediction quality
is visible. Informational (not a pass/fail gate): the union spans heterogeneous
platforms and ~13x more samples, so absolute metric levels are not comparable.

Usage:
    python 05_xpert_finetune_bench/compare_v1_benchmark.py
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, V1_WORKFLOW_ROOT  # noqa: E402

METRICS = ["MSE_deg", "RMSE_deg", "MAE_deg",
           "Pearson_deg", "Spearman_deg", "R2_deg"]
OUT = WORKFLOW_ROOT / "05_xpert_finetune_bench" / "output" / "benchmark_vs_v1.csv"


def load_bench(root, bench_step, split="warm-start"):
    p = root / bench_step / "output" / "benchmark.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    return df[df["split"] == split]


def main():
    v2 = load_bench(WORKFLOW_ROOT, "05_xpert_finetune_bench")
    # v1's benchmark step is 04_xpert_finetune_bench (pre-union layout).
    v1 = load_bench(V1_WORKFLOW_ROOT, "04_xpert_finetune_bench")
    if v2 is None:
        print("ERROR: v2 benchmark.csv not found — run xpert_finetune_bench.py first")
        sys.exit(1)

    rows = []
    for m in METRICS:
        rows.append({
            "metric": m,
            "v2_union_mean": v2[m].mean() if m in v2 else None,
            "v2_union_std": v2[m].std() if m in v2 else None,
            "v1_dtoxs_mean": v1[m].mean() if v1 is not None and m in v1 else None,
            "v1_dtoxs_std": v1[m].std() if v1 is not None and m in v1 else None,
        })
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(out.to_string(index=False))
    print(f"\nSaved: {OUT}")
    if v1 is None:
        print("(v1 benchmark.csv not found — v1 columns left empty)")


if __name__ == "__main__":
    main()
