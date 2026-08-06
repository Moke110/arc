#!/usr/bin/env python3
"""
Step 01 — Drug list.

Generates the 3,127-drug FDA list consumed by every downstream step. The
source of truth is the ChEMBL-derived FDA drug table (input/fda_drugs.csv).
Output is written to 01_drug_list/output/drug_list.csv with the columns:
chembl_id, drug_name, canonical_smiles, drug_idx. drug_idx is the shared
11,104-drug embedding index (see input/all_drugs_idx2smi.npy); it is resolved
from the precomputed SMILES mapping when available.

Usage:
    python 01_drug_list/01_drug_list.py
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "01_drug_list"


def resolve_drug_indices(drug_df, smi_path):
    """Attach the shared embedding drug_idx to each drug via its SMILES."""
    idx2smi = np.load(smi_path, allow_pickle=True).item()
    smi2idx = {smi: idx for idx, smi in idx2smi.items()}
    drug_df = drug_df.copy()
    drug_df["drug_idx"] = drug_df["canonical_smiles"].map(smi2idx)
    n_missing = drug_df["drug_idx"].isna().sum()
    if n_missing:
        print(f"WARNING: {n_missing} drugs without a precomputed embedding index "
              f"(drug_idx will be -1)")
        drug_df["drug_idx"] = drug_df["drug_idx"].fillna(-1).astype(int)
    return drug_df


def main():
    parser = argparse.ArgumentParser(description="Generate the FDA drug list")
    parser.add_argument("--drug-csv",
                        default=str(WORKFLOW_ROOT / "input/fda_drugs.csv"),
                        help="Source FDA drug table (chembl_id, pref_name, canonical_smiles)")
    parser.add_argument("--output",
                        default=str(STEP_DIR / "output/drug_list.csv"),
                        help="Output drug list CSV")
    parser.add_argument("--smi-mapping",
                        default=str(WORKFLOW_ROOT / "input/all_drugs_idx2smi.npy"),
                        help="Precomputed idx->SMILES mapping (.npy)")
    args = parser.parse_args()

    drug_df = pd.read_csv(args.drug_csv)
    drug_df = drug_df.rename(columns={"pref_name": "drug_name"})
    if "drug_name" not in drug_df.columns:
        drug_df["drug_name"] = drug_df.iloc[:, 1]
    drug_df = resolve_drug_indices(drug_df, args.smi_mapping)

    out_cols = [c for c in ["chembl_id", "drug_name", "canonical_smiles", "drug_idx"]
                if c in drug_df.columns]
    out = drug_df[out_cols]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"Step 01 complete — {len(out)} drugs -> {args.output}")
    print(f"  Drugs with embedding index: {(out['drug_idx'] >= 0).sum()}/{len(out)}")


if __name__ == "__main__":
    main()
