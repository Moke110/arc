"""
Shared paths and constants for the Stage 2 v1 CM benchmark workflow
(`stage2/_archive/cm_bench_v1/`).

Stage 2 benchmarks **whole-transcriptome** drug-perturbation models on
cardiomyocytes (ADR 0020). All paths resolve from the workflow root so every
numbered stage can run from any working directory. Raw CM datasets under
`pert_data/` are read-only; every artifact is written under this workflow.
"""

from pathlib import Path

# ── Roots ─────────────────────────────────────────────────────────────
WORKFLOW_ROOT = Path(__file__).resolve().parent.parent
WORKSPACE_ROOT = WORKFLOW_ROOT.parent.parent          # repository root
CM_DATA_ROOT = WORKSPACE_ROOT / "pert_data"
INPUT_DIR = WORKFLOW_ROOT / "input"
OUTPUT_DIR = WORKFLOW_ROOT / "output"
BENCH_DIR = OUTPUT_DIR / "bench"

# ── Benchmark scope ───────────────────────────────────────────────────
# The three participating whole-transcriptome CM datasets. gse262419
# (TempO-Seq, mostly environmental compounds) is excluded by the spec.
DATASETS = ["dtoxs_primary", "dtoxs_ipsc", "gse244740"]

# Identifier space of each dataset's raw gene index.
GENE_ID_SPACE = {
    "dtoxs_primary": "symbol",
    "dtoxs_ipsc": "symbol",
    "gse244740": "ensembl",
}

# Controls per dataset (dataset-wide control states / plate controls).
DTOXS_CONTROL_STATES = {"CTRL", "CTRL-DMSO", "CTRL-EtOH"}

# Dose encoding is log10(uM), shared across datasets (spec: Model I/O contract).
DOSE_EPS = 1e-3  # uM floor so log10 is finite for the smallest doses

# A (dataset, drug, dose) group is eligible for scoring only with at least this
# many replicate treated wells (spec: Ground truth).
MIN_WELLS_PER_GROUP = 4

# Ground-truth significance threshold.
DEG_FDR = 0.05
DEG_ABS_MEAN_GATE = 0.0  # optional |mean Δ| gate; 0 disables it

# Cold-drug CV folds.
COLD_DRUG_FOLDS = 5
COLD_DRUG_SEED = 0

# cold-dose axis (gse244740 only): one interior + one extreme level held out.
COLD_DOSE_INTERIOR_LEVEL = 4   # Con4 (~ middle of the 8-point series)
COLD_DOSE_EXTREME_LEVEL = 8    # Con8 (top concentration)


def dataset_path(dataset):
    return CM_DATA_ROOT / dataset
