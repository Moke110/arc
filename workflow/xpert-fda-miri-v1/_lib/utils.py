"""
Shared utilities for the finalized xpert-fda-miri-v1 MIRI screening workflow.

xpert-fda-miri-v1 is the **author-aligned, L1000-contained FDA screen**: it
predicts only the FDA drugs already present in the L1000 embedding table
(exact canonical-SMILES match, drug_idx < 8981), so no HG training /
fine-tuning / FDA-extension inference is ever needed. The HG biological
embeddings are the author reference table itself, workflow-owned in
input/HG_drug_embeddings.npy. All paths in this module are resolved from the
workflow root so each numbered stage can be run independently from any working
directory.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Workflow / component roots ────────────────────────────────────────
WORKFLOW_ROOT = Path(__file__).resolve().parent.parent
WORKSPACE_ROOT = WORKFLOW_ROOT.parent.parent          # repository root
XPERT_ROOT = WORKSPACE_ROOT / "pert_model" / "xpert"
DISEASE_DATA_ROOT = WORKSPACE_ROOT / "disease_data" / "miri"
# The v1-era DToxS-primary well-level package, frozen byte-identically at
# pert_data/dtoxs_primary/v1-frozen/ (the data xpert-fda-miri-v1 used; the
# live LINCS LDS-1593 package now lives flat at pert_data/dtoxs_primary/).
DTOXS_ROOT = WORKSPACE_ROOT / "pert_data" / "dtoxs_primary" / "v1-frozen"

# L1000 embedding-table boundary: rows [0, N_L1000_DRUGS) are the author
# L1000 reference; rows >= N_L1000_DRUGS are the FDA-extension rows and are
# OUT OF SCOPE for xpert-fda-miri-v1 (never touched).
N_L1000_DRUGS = 8981

# Author HG reference, workflow-owned in input/ (moved here from the HG
# component so the workflow is self-contained on drug assets).
AUTHOR_HG_REFERENCE = WORKFLOW_ROOT / "input" / "HG_drug_embeddings.npy"

# Insert the XPert model library so `python -m src.*` and `from src.*`
# resolve (HG is a pure library and is never imported by this workflow).
sys.path.insert(0, str(XPERT_ROOT))


def ensure_arc_env():
    """Best-effort check that the arc conda env provides the core stack."""
    try:
        import torch  # noqa: F401
    except ImportError:
        print("WARNING: arc conda env not active — activate it first "
              "(conda activate arc)")


def load_drug_deg(deg_csv_path, zscore=False):
    """Load a predicted drug DEG matrix (drug_name, smiles, drug_idx + 978 genes).

    Returns (meta DataFrame, gene_cols list, deg_matrix ndarray).
    """
    df = pd.read_csv(deg_csv_path)
    gene_cols = df.columns[3:].tolist()  # skip drug_name, smiles, drug_idx
    meta = df[["drug_name", "smiles", "drug_idx"]].copy()
    deg_matrix = df[gene_cols].values.astype(float)
    if zscore:
        deg_matrix = zscore_normalize(deg_matrix)
    return meta, gene_cols, deg_matrix


def zscore_normalize(deg_matrix):
    """Z-score normalize per gene across drugs."""
    gene_mean = deg_matrix.mean(axis=0)
    gene_std = deg_matrix.std(axis=0)
    gene_std[gene_std < 1e-10] = 1.0
    return (deg_matrix - gene_mean) / gene_std


def resolve_l1000_drug_rows(drug_df):
    """Filter a drug table to the L1000-contained rows (drug_idx < 8981).

    Returns (kept_df, dropped_df) where dropped holds drugs whose exact
    canonical SMILES has no match in the L1000 embedding table.
    """
    kept = drug_df[drug_df["drug_idx"].between(0, N_L1000_DRUGS - 1)].copy()
    dropped = drug_df[~drug_df["drug_idx"].between(0, N_L1000_DRUGS - 1)].copy()
    return kept, dropped
