#!/usr/bin/env python3
"""
Stage 08 — Efficacy scores.

Runs the three MIRI efficacy scorers (WCS+RGES, Zhang/ssCMap, pathway RGES)
against the integrated MIRI reference and the curated pathway gene sets, and
computes the four z-components + integrated score for all drugs.

Output: data/scores/1-scores.csv

Usage:
    python scripts/08_efficacy_scores.py
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT  # noqa: E402
from scoring import run_scoring  # noqa: E402

DEG = WORKFLOW_ROOT / "data" / "deg-zs" / "xpert_fda-cm_deg-zs.csv"
GENE_SET = WORKFLOW_ROOT / "reference-genesets" / "MIRI_up-down.json"
RANK_REF = WORKFLOW_ROOT / "reference-genesets" / "MIRI_rank.json"
PATHWAYS = WORKFLOW_ROOT / "reference-genesets" / "MIRI_pathways.json"
OUTPUT = WORKFLOW_ROOT / "data" / "scores" / "1-scores.csv"


def main():
    parser = argparse.ArgumentParser(description="Compute MIRI efficacy scores")
    parser.add_argument("--deg-csv", default=str(DEG))
    parser.add_argument("--gene-set", default=str(GENE_SET))
    parser.add_argument("--rank-ref", default=str(RANK_REF))
    parser.add_argument("--pathways", default=str(PATHWAYS))
    parser.add_argument("--output", default=str(OUTPUT))
    parser.add_argument("--no-zscore", action="store_true",
                        help="Disable per-gene z-score normalization")
    parser.add_argument("--deg-zs", action="store_true", default=True,
                        help="Input DEG is already z-scored (stage 06 output); "
                             "do not re-normalize. Pass --no-deg-zs to score the "
                             "raw predicted DEG (z-scoring applied inside).")
    args = parser.parse_args()

    zscore = not args.no_zscore and not args.deg_zs
    scores = run_scoring(args.deg_csv, args.gene_set, args.rank_ref,
                         args.pathways, zscore=zscore)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.output, index=False)

    print(f"Stage 08 complete — {len(scores)} drugs scored")
    print(f"  Saved to {args.output}")
    print(f"\n  Top 10 by integrated score:")
    for i, (_, row) in enumerate(scores.head(10).iterrows()):
        print(f"    {i+1:2d}. {row['drug_name']:40s} integrated={row['integrated']:+.4f}")


if __name__ == "__main__":
    main()
