"""
Shared utilities for the xpert-fda-miri-v2 MIRI screening workflow.

xpert-fda-miri-v2 is the **union-CM finetune sibling** of the finalized
author-aligned L1000-contained FDA screen (xpert-fda-miri-v1). It screens the
same 1,004 FDA drugs already present in the L1000 embedding table (exact
canonical-SMILES match, drug_idx < 8981) with the same baseline, MIRI scoring,
filters, integration and report — the ONLY substantive change is the XPert
finetune data: instead of the frozen DToxS-only input/dtoxs.h5ad, steps
05/06/07 train on the workflow-built **union h5ad** of all five CM datasets
under pert_data/cm (built by 03_cm_finetune_data). The HG biological
embeddings are the author reference table, workflow-owned in
input/HG_drug_embeddings.npy.

All paths in this module are resolved from the workflow root so each numbered
stage can be run independently from any working directory.

Note on the v1 stale constant: v1's utils resolved DTOXS_ROOT =
pert_data/dtoxs, a path that no longer exists after the 2026-08-12 rename
(dtoxs -> dtoxs-xpert + dtoxs_primary). v2 has no DTOXS_ROOT; the union prep
step reads per-dataset paths directly under pert_data/cm, and the step-02
baseline reads the raw GEO deposit under pert_data/cm/dtoxs_primary.
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

# CM finetune-data union sources (read-only; never modified by the workflow).
CM_DATA_ROOT = WORKSPACE_ROOT / "pert_data" / "cm"
CM_DATASETS = ["dtoxs_primary", "dtoxs_ipsc", "gse244740", "gse247541",
               "gse262419"]

# The finalized v1 workflow (read-only schema reference for the union gate and
# the DToxS baseline comparison). Only ever read, never written.
V1_WORKFLOW_ROOT = WORKSPACE_ROOT / "workflow" / "xpert-fda-miri-v1"
V1_DTOXS_H5AD = V1_WORKFLOW_ROOT / "input" / "dtoxs.h5ad"

# L1000 embedding-table boundary: rows [0, N_L1000_DRUGS) are the author
# L1000 reference; rows >= N_L1000_DRUGS are the FDA-extension rows and are
# OUT OF SCOPE for xpert-fda-miri-v1/v2 (never touched).
N_L1000_DRUGS = 8981

# Author HG reference, workflow-owned in input/ (moved here from the HG
# component so the workflow is self-contained on drug assets).
AUTHOR_HG_REFERENCE = WORKFLOW_ROOT / "input" / "HG_drug_embeddings.npy"

# Author exposure-time map (the XPert-author L1000 convention; identical to
# the pretrain and to v1's frozen dtoxs.h5ad, which is all 24 h -> index 2).
# The union finetune data is restricted to these times (step 03), so the
# 48 h / 144 h datasets (dtoxs_primary, dtoxs_ipsc, gse247541) are excluded.
AUTHOR_TIME_MAP = {3: 0, 6: 1, 24: 2}
AUTHOR_TIMES = (3, 6, 24)

# L1000 XPert pretrain h5ad — the pooled log2(CPM+1) distribution used as the
# quantile-alignment reference for the union finetune data (step 03).
L1000_PRETRAIN_H5AD = (CM_DATA_ROOT.parent / "l1000-xpert" / "data"
                       / "l1000_mdmt_268022_pretrain.h5ad")

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
