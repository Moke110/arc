"""
Shared utilities for the xpert-fda-v1 MIRI repositioning workflow.

The workflow owns every script that generates its artifacts; the model
components (pert_model/xpert, drug_emb_model/hg) are imported as libraries.
All paths in this module are resolved from the workflow root so each numbered
stage can be run independently from any working directory.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Workflow / component roots ────────────────────────────────────────
WORKFLOW_ROOT = Path(__file__).resolve().parent.parent
WORKSPACE_ROOT = WORKFLOW_ROOT.parent.parent          # repository root
XPERT_ROOT = WORKSPACE_ROOT / "pert_model" / "xpert"
HG_ROOT = WORKSPACE_ROOT / "drug_emb_model" / "hg"
DISEASE_DATA_ROOT = WORKSPACE_ROOT / "disease_data" / "miri"

# Insert the model library so `python -m src.*` and `from src.*` resolve.
sys.path.insert(0, str(XPERT_ROOT))
sys.path.insert(0, str(HG_ROOT))


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
