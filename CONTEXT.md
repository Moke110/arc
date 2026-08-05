# ARC

AI-driven drug repositioning for cardiac diseases. Predict the transcriptional perturbation a drug causes on cardiomyocytes, then rank drugs by predicted MIRI efficacy. The repository is split into three categories: dev-tool files, drug perturbation prediction components, and drug repositioning workflows.

## Language

### Categories

**Dev-tool files**:
Files that support the workspace itself rather than the science: `.gitignore`, `skills-lock.json`, `.agents/`, `docs/`, `reference/`, `AGENTS.md`, `environment.yml`.
_Avoid_: tooling, meta

**Perturbation component**:
One of the three building blocks of prediction: a perturbation dataset (`pert_data/`), a drug embedding model (`drug_emb_model/`), or a perturbation prediction model (`pert_model/`).
_Avoid_: module, stage

**Repositioning workflow**:
A complete end-to-end drug-screening pipeline for one disease, bundling its scripts, prompts, data, and notebook. `miri_reposition/xpert-fda-v1/` is the MIRI workflow.
_Avoid_: pipeline (when meaning a workflow), project

**Disease data**:
Shared disease-specific inputs, reusable across future workflows: raw GEO acquisitions and per-GSE DEG lists under `disease_data/`. Distinct from workflow artifacts because more than one workflow can consume it.
_Avoid_: input data, shared data

### Components

**Drug embedding**:
A learned vector representation of a drug: HG biological embedding (256-dim) and UniMol 3D features (122×514). Produced by embedding-model inference; consumed by finetune and predict.
_Avoid_: drug feature

**Perturbed transcriptome**:
The predicted expression of 978 L1000 landmark genes for a drug-treated cell. Produced by XPert inference; z-scored and scored in the workflow.

**Drug list**:
The set of drugs a workflow screens — for `xpert-fda-v1`, 3,127 FDA-approved drugs. Identity of each drug is its SMILES; the canonical mapping is `all_drugs_idx2smi.npy`.
_Avoid_: drug set

**Baseline transcriptome**:
The consensus control expression (978 genes) of the cell type under study, derived from DTOXS control wells. Workflow-owned; the reference condition the perturbed transcriptome is measured against.

### Workflow stages

**MIRI efficacy score**:
A score of how strongly a drug's predicted perturbation reverses the MIRI transcriptional disease signature. Computed by three methods (WCS+RGES, Zhang/ssCMap, pathway RGES) yielding four z-components.
_Avoid_: drug score

**Cardiotoxicity score**:
Deferred from v1. A score of a drug's risk of cardiac harm (e.g., hERG/QT).
_Avoid_: safety score

**MIRI reference DEG list**:
The integrated, cross-dataset MIRI differential-expression reference (`MIRI_up-down.json`, `MIRI_rank.json`) built from the per-GSE DEG lists. Distinct from per-GSE DEG lists and from the curated `MIRI_pathways.json`.

**Per-GSE DEG list**:
The DEG list derived from a single GEO study, kept alongside its acquisition data in `disease_data/miri/{GSE}/deg/`. The raw material for the MIRI reference DEG list.

### Boundary

**Producer-as-library**:
Components ship library code (`src/`) and trained weights; workflows hold the scripts that import the libraries and generate every workflow artifact. A component does not write into a workflow's folder at runtime; workflows regenerate artifacts by running their scripts.
_Avoid_: dependency, integration
