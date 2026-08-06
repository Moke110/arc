#!/usr/bin/env python3
"""
Toolchain-parity verification against v1 (gate E).

Because v2's rows 0-8980 are byte-identical to the author reference (gate B),
the unchanged L1000 side must produce unchanged training-side results. Asserts:

  1. Steps 01 and 03 regenerate v1's drug list and baseline transcriptome
     byte-identically (the workflow root resolves to v2 but the computation is
     the same).
  2. Step 08 regenerates v1's MIRI reference gene sets byte-identically
     (workflow-independent).
  3. The v2 DToxS benchmark (`05_benchmark_xpert/output/benchmark.csv`)
     matches v1's per-fold metrics within float precision — Pearson and
     Spearman must match exactly; MSE/RMSE/MAE within 1e-2 (GPU-fp
     nondeterminism in the finetune, documented in AGENTS.md).

Downstream gate (observational per the feature spec): the FDA screening tail
(steps 06-11) is expected to diverge from v1 — captured by
`11_integrate/divergence_vs_v1.md` (`_lib/compare_v1.py`), not asserted here.

Usage:
    python _lib/verify_parity_v1.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, WORKSPACE_ROOT  # noqa: E402

V1 = WORKSPACE_ROOT / "miri_reposition" / "xpert-fda-v1"

CSV_PAIRS = [
    ("01_drug_list/output/drug_list.csv", "01_drug_list/output/drug_list.csv"),
    ("03_baseline_transcriptome/output/cm_ctl_landmark.csv",
     "03_baseline_transcriptome/output/cm_ctl_landmark.csv"),
    ("08_miri_ref_path/output/MIRI_up-down.json",
     "08_miri_ref_path/output/MIRI_up-down.json"),
    ("08_miri_ref_path/output/MIRI_rank.json",
     "08_miri_ref_path/output/MIRI_rank.json"),
]

results = []


def check(name, cond, detail=""):
    ok = bool(cond)
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def verify_unchanged_outputs():
    print("\n1. Unchanged steps reproduce v1 outputs")
    ok = True
    for rel_v1, rel_v2 in CSV_PAIRS:
        p1, p2 = V1 / rel_v1, WORKFLOW_ROOT / rel_v2
        if not p1.exists() or not p2.exists():
            ok &= check(f"{rel_v2} both exist", False,
                        f"{p1.exists()} / {p2.exists()}")
            continue
        if p1.suffix == ".json":
            same = p1.read_text() == p2.read_text()
        else:
            same = pd.read_csv(p1).equals(pd.read_csv(p2))
        ok &= check(f"{rel_v2} identical to v1", same)
    return ok


def verify_benchmark_parity():
    print("\n2. DToxS benchmark matches v1 (toolchain parity on unchanged L1000 rows)")
    v1 = pd.read_csv(V1 / "05_benchmark_xpert/output/benchmark.csv")
    v2 = pd.read_csv(WORKFLOW_ROOT / "05_benchmark_xpert/output/benchmark.csv")
    if not (v1[["split", "fold"]].values == v2[["split", "fold"]].values).all():
        return check("same split/fold rows", False)

    ok = True
    for col in ["Pearson", "Spearman"]:
        ok &= check(f"{col} identical", np.array_equal(v1[col], v2[col]))
    for col in ["MSE", "RMSE", "MAE", "MSE_deg", "RMSE_deg", "MAE_deg"]:
        d = np.abs(v1[col].values.astype(float) - v2[col].values.astype(float)).max()
        ok &= check(f"{col} within float precision", d < 1e-2, f"max diff {d:.2e}")
    return ok


def main():
    print("=" * 60)
    print("v2 toolchain parity vs v1 (gate E)")
    print("=" * 60)
    ok = verify_unchanged_outputs()
    ok = verify_benchmark_parity() and ok
    print("\n" + "=" * 60)
    print(f"RESULT: {'PARITY CONFIRMED' if ok else 'PARITY BROKEN'}")
    print("=" * 60)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
