#!/usr/bin/env python3
"""
Step 01 — Drug list (xpert-fda-miri-v1: L1000-contained FDA drugs only).

The finalized MIRI screen predicts only FDA drugs that are already present in
the L1000 embedding table — i.e. drugs whose **exact canonical SMILES** appears
in the L1000 rows (drug_idx in [0, 8981)) of the shared idx->SMILES mapping
(input/all_drugs_idx2smi.npy). Drugs matching only the FDA-extension rows
(8981-11103) or unmatched are dropped: the screen never needs HG training /
FDA-extension inference because every screened drug uses the author-reference
HG embedding directly.

Output: 01_drug_embed/output/drug_list.csv
  Columns: chembl_id, drug_name, canonical_smiles, drug_idx
  Rows: the L1000-contained FDA drugs (exact-SMILES match, drug_idx < 8981).

Usage:
    python 01_drug_embed/drug_list.py
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, N_L1000_DRUGS, resolve_l1000_drug_rows  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "01_drug_embed"


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
    parser = argparse.ArgumentParser(description="Generate the L1000-contained FDA drug list")
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

    kept, dropped = resolve_l1000_drug_rows(drug_df)
    print(f"FDA drugs total: {len(drug_df)}")
    print(f"  in L1000 (drug_idx < {N_L1000_DRUGS}, exact SMILES match): {len(kept)}")
    print(f"  FDA-extension rows / unmatched (dropped for xpert-fda-miri-v1): {len(dropped)}")

    out_cols = [c for c in ["chembl_id", "drug_name", "canonical_smiles", "drug_idx"]
                if c in kept.columns]
    kept = kept[out_cols]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    kept.to_csv(args.output, index=False)
    print(f"Step 01 complete — {len(kept)} drugs -> {args.output}")


if __name__ == "__main__":
    main()
