#!/usr/bin/env python3
"""
Stage 10 — Filter + rank.

Applies the screening filters and produces the final ranked drug list:

  1. WCS filter:      keep drugs with WCS < 0
  2. Pathway filter:  keep drugs with >= MIN_CRITICAL critical pathways whose
                      efficacy > 0
  3. Rank:            sort the survivors by integrated score descending

The cardiotoxicity stage is deferred from v1 (see docs/adr/0004); no
cardiotoxicity filter is applied.

Output:
  data/ranked/ranked_drugs.csv        — filtered + ranked drugs (final result)
  data/ranked/2-wcs_filtered.csv      — WCS < 0 drugs
  data/ranked/3-path_filtered.csv     — pathway-filtered drugs

Usage:
    python scripts/10_filter_rank.py
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT  # noqa: E402
from scoring import load_critical_names  # noqa: E402

SCORES = WORKFLOW_ROOT / "data" / "scores" / "1-scores.csv"
PATHWAYS = WORKFLOW_ROOT / "reference-genesets" / "MIRI_pathways.json"
RANKED_OUT = WORKFLOW_ROOT / "data" / "ranked" / "ranked_drugs.csv"
MIN_CRITICAL = 6


def main():
    parser = argparse.ArgumentParser(description="Filter + rank MIRI drug candidates")
    parser.add_argument("--scores", default=str(SCORES))
    parser.add_argument("--pathways", default=str(PATHWAYS))
    parser.add_argument("--output", default=str(RANKED_OUT))
    parser.add_argument("--min-critical", type=int, default=MIN_CRITICAL,
                        help="Minimum critical pathways with efficacy > 0")
    args = parser.parse_args()

    scores = pd.read_csv(args.scores)
    critical = load_critical_names(args.pathways)
    print(f"Critical pathways ({len(critical)}):")

    # Step 1: WCS < 0
    wcs_pass = scores[scores["WCS"] < 0].copy()
    print(f"  WCS < 0: {len(wcs_pass)}/{len(scores)} drugs")

    # Step 2: pathway filter
    n_positive = pd.Series(0, index=wcs_pass.index)
    for name in critical:
        col = f"{name}_efficacy"
        if col not in wcs_pass.columns:
            print(f"  WARNING: missing column {col} — skipping pathway filter")
            continue
        n_positive += (wcs_pass[col] > 0).astype(int)
    path_pass = wcs_pass[n_positive >= args.min_critical].copy()
    print(f"  >= {args.min_critical}/{len(critical)} critical pathways > 0: "
          f"{len(path_pass)}/{len(wcs_pass)} drugs")

    # Step 3: rank by integrated descending
    ranked = path_pass.sort_values("integrated", ascending=False).reset_index(drop=True)

    out_dir = Path(args.output).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    Path(out_dir / "2-wcs_filtered.csv").write_text(
        wcs_pass[["drug_name", "integrated"]].sort_values(
            "integrated", ascending=False).to_csv(index=False))
    Path(out_dir / "3-path_filtered.csv").write_text(
        ranked[["drug_name", "integrated"]].to_csv(index=False))
    ranked.to_csv(args.output, index=False)

    print(f"\nStage 10 complete — {len(ranked)} drugs ranked")
    print(f"  Saved: {args.output}")
    print(f"\n  Top candidates:")
    for i, (_, row) in enumerate(ranked.head(10).iterrows()):
        print(f"    {i+1:2d}. {row['drug_name']:40s} integrated={row['integrated']:+.4f}")


if __name__ == "__main__":
    main()
