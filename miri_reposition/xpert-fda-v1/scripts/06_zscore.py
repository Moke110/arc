#!/usr/bin/env python3
"""
Stage 06 — Z-score of perturbed transcriptome.

Z-score normalizes the predicted DEG matrix per gene across drugs, removing the
dominant gene-level baseline so inter-drug differences drive downstream
scoring. Always run with this stage (AGENTS.md: ranking scripts take --zscore).

Output: data/deg-zs/xpert_fda-cm_deg-zs.csv

Usage:
    python scripts/06_zscore.py
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, zscore_normalize  # noqa: E402

INPUT = WORKFLOW_ROOT / "data" / "deg" / "xpert_fda-cm_deg.csv"
OUTPUT = WORKFLOW_ROOT / "data" / "deg-zs" / "xpert_fda-cm_deg-zs.csv"


def main():
    parser = argparse.ArgumentParser(description="Z-score normalize the predicted DEG matrix")
    parser.add_argument("--deg-csv", default=str(INPUT), help="Input predicted DEG CSV")
    parser.add_argument("--output", default=str(OUTPUT), help="Output z-scored DEG CSV")
    args = parser.parse_args()

    df = pd.read_csv(args.deg_csv)
    gene_cols = df.columns[3:].tolist()
    z = zscore_normalize(df[gene_cols].values.astype(float))

    out = df[["drug_name", "smiles", "drug_idx"]].copy()
    out[gene_cols] = z
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    print(f"Stage 06 complete — {len(out)} drugs x {len(gene_cols)} genes z-scored")
    print(f"  Saved to {args.output}")


if __name__ == "__main__":
    main()
