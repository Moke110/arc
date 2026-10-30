"""
Shared utilities for the xpert-fda-miri-v3 MIRI screening workflow.

xpert-fda-miri-v3 is the **all-CM finetune sibling** of the finalized
author-aligned L1000-contained FDA screen (xpert-fda-miri-v1, ADR 0007). It
screens the same 1,004 FDA drugs already present in the L1000 embedding table
(exact canonical-SMILES match, drug_idx < 8981) with the same baseline, MIRI
scoring, filters and integration — the ONLY substantive change is the XPert
finetune data: instead of v1's frozen DToxS-only finetune h5ad
(`stage1/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad`), steps
05/06/07 train on the workflow-built **union h5ad of all four usable CM
datasets** (dtoxs_primary + dtoxs_ipsc + gse244740 + gse262419; ADR 0017),
per-plate normalized (ADR 0016). gse247541 is excluded (6-day exposure outside
the author time map, 2 in-HG drugs). v2 stays authoritative and untouched.

All paths in this module are resolved from the workflow root so each numbered
stage can be run independently from any working directory.

Data-lake layout (flattened 2026-08-17): the CM datasets live directly under
pert_data/ (the `cm/` grouping was removed), i.e. pert_data/LINCS-LDS-1593/source,
pert_data/DToxS-iPSC/source, pert_data/GSE244740/source, pert_data/GSE247541/source,
pert_data/GSE262419/source, plus pert_data/XPert-L1000/source (pretrain corpus). The DToxS
primary study has one deposit under pert_data/LINCS-LDS-1593/source: the released
LINCS LDS-1593 well-level package (Data/, Metadata/, ...; the GSE146096/146095
raw GEO deposit was deleted 2026-08-18 as a strict subset), which carries the
non-KI perturbations and the DToxS small-molecule metadata. The v1-era
well-level package is frozen byte-identically at
pert_data/LINCS-LDS-1593/local/v1-frozen/. The finetune-data prep reads the LINCS
wells; it never modifies a source file.
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
CM_DATA_ROOT = WORKSPACE_ROOT / "pert_data"
CM_DATASETS = ["dtoxs_primary", "dtoxs_ipsc", "gse244740", "gse247541",
               "gse262419"]
# The step-03 finetune scope since ADR 0017: all four usable CM datasets
# (primary + iPSC + gse244740 + gse262419), per-plate normalized (ADR 0016).
# gse247541 is excluded (6-day exposure outside the author time map, 2 in-HG
# drugs, 12 treated samples) but stays in the data lake + drug-overlap audit.
FINETUNE_DATASETS = ["dtoxs_primary", "dtoxs_ipsc", "gse244740", "gse262419"]

# The released LINCS LDS-1593 well-level package (DToxS primary, Level 3),
# now living flat at pert_data/LINCS-LDS-1593/source/: per-batch well-level
# read-counts + experiment designs + the LINCS metadata suite
# (Small_Molecule_Metadata.txt etc.). Released 2026-08-18 (the GSE146095/
# GSE146096 raw GEO deposit was deleted; the package is its strict superset).
# The v1-era well-level package is frozen byte-identically at
# pert_data/LINCS-LDS-1593/local/v1-frozen/ (the data xpert-fda-miri-v1 used).
DTOXS_LINCS_ROOT = CM_DATA_ROOT / "LINCS-LDS-1593" / "source"
DTOXS_LINCS_LEVEL3 = DTOXS_LINCS_ROOT / "Data" / "Level3"
DTOXS_SMALL_MOL = DTOXS_LINCS_ROOT / "Metadata" / "Small_Molecule_Metadata.txt"

# The finalized v1 workflow (read-only schema reference for the union gate and
# the DToxS baseline comparison). Only ever read, never written.
V1_WORKFLOW_ROOT = WORKSPACE_ROOT / "workflow" / "xpert-fda-miri-v1"
V1_DTOXS_H5AD = V1_WORKFLOW_ROOT / "03_cm_finetune_data" / "output" / "dtoxs.h5ad"

# The DToxS-LINCS testbed (v1.1): its rebuilt finetune h5ad is the benchmark
# reference the v3 primary arm realigns to (ADR 0015). v3 imports v1.1's fold
# columns verbatim for primary rows (byte-identical test sets, ADR 0017). The
# v2 union (per-plate, ADR 0016) is the delta/composition reference the four-
# dataset union is audited against. Only ever read.
V11_WORKFLOW_ROOT = WORKSPACE_ROOT / "workflow" / "xpert-fda-miri-v1.1"
V11_DTOXS_LINCS_H5AD = (V11_WORKFLOW_ROOT / "03_cm_finetune_data" / "output"
                        / "dtoxs_lincs.h5ad")
V2_WORKFLOW_ROOT = WORKSPACE_ROOT / "workflow" / "xpert-fda-miri-v2"
V2_CM_UNION_H5AD = (V2_WORKFLOW_ROOT / "03_cm_finetune_data" / "output"
                    / "cm_union.h5ad")

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
# quantile-alignment reference for the union finetune data (step 03). The data
# lake lays the pretrain corpus out under the XPert-L1000 item's source/ tree.


def l1000_data_h5ad(filename):
    """Resolve one L1000 corpus h5ad under the XPert-L1000 item's source/ tree."""
    return CM_DATA_ROOT / "XPert-L1000" / "source" / "data" / filename


L1000_PRETRAIN_H5AD = l1000_data_h5ad("l1000_mdmt_268022_pretrain.h5ad")
L1000_FULL_H5AD = l1000_data_h5ad("l1000_mdmt_full_336852.h5ad")

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
