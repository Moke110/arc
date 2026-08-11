# ARC

AI-driven drug repositioning for cardiac diseases. Predict the transcriptional perturbation a drug causes on cardiomyocytes, then rank drugs by predicted MIRI efficacy. The repository is split into three categories: dev-tool files, drug perturbation prediction components, and drug repositioning workflows.

## Language

### Categories

**Dev-tool files**:
Files that support the workspace itself rather than the science: `.gitignore`, `skills-lock.json`, `.agents/`, `docs/`, `reference/`, `AGENTS.md`, `environment.yml`.
_Avoid_: tooling, meta

**Perturbation component**:
One of the three building blocks of prediction: a perturbation dataset (`pert_data/`), a drug embedding model (`drug_emb_model/`), or a perturbation prediction model (`pert_model/`). Components are pure static libraries (ADR 0008): `src/` code + small metadata only; data/weights live in the consuming workflow or dataset home.
_Avoid_: module, stage

**Repositioning workflow**:
A complete end-to-end drug-screening pipeline for one disease, bundling its scripts, prompts, data, and notebook. MIRI workflows live under `workflow/`: `workflow/xpert-fda-miri-v1/` is the canonical finalized MIRI screen; `workflow/xpert-fda-miri-v0/` is the frozen legacy snapshot.
_Avoid_: pipeline (when meaning a workflow), project

**Disease data**:
Shared disease-specific inputs, reusable across future workflows: raw GEO acquisitions and per-GSE DEG lists under `disease_data/`. Distinct from workflow artifacts because more than one workflow can consume it.
_Avoid_: input data, shared data

### Components

**Drug embedding**:
A learned vector representation of a drug: HG biological embedding (256-dim) and UniMol 3D features (122×514). Produced by embedding-model inference; consumed by finetune and predict. The finalized MIRI screen uses only the author-reference HG table (byte-identical) plus the precomputed UniMol array.
_Avoid_: drug feature

**Perturbed transcriptome**:
The predicted expression of 978 L1000 landmark genes for a drug-treated cell. Produced by XPert inference; z-scored and scored in the workflow.

**Drug list**:
The set of drugs a workflow screens — for `xpert-fda-miri-v1`, the 1,004 L1000-contained FDA drugs (exact canonical-SMILES match, `drug_idx < 8981`). Identity of each drug is its SMILES; the canonical mapping is `all_drugs_idx2smi.npy`.
_Avoid_: drug set

**Baseline transcriptome**:
The consensus control expression (978 genes) of the cell type under study, derived from DTOXS control wells. Workflow-owned; the reference condition the perturbed transcriptome is measured against.

### Workflow stages (finalized step layout)

The MIRI screen reads top-down: `01_drug_embed`, `02_cm_baseline`, `03_xpert_pretrain`, `04_xpert_finetune_bench`, `05_xpert_finetune`, `06_xpert_predict`, `07_miri_phenotype`, `08_miri_efficacy_scores`, `09_filter`, `10_integrate`, `11_report`.

**Pretrain unit**:
The workflow-owned author pretrained XPert weight + pretrain YAML (`03_xpert_pretrain/`); the L1000 pretrain h5ads live in `pert_data/l1000/` (dataset home). Both pretrain and finetune YAML recipes are owned by the workflow.

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

**Final drug report**:
The step-11 agentic bilingual (EN/CN) per-drug fact-sheet report over the `09_filter` survivors with integrated-score percentile > 0.5 (`11_report/`), driven by `prompt.md` and assembled from resumable per-drug agent JSONs.
_Avoid_: drug research report

### Boundary

**Producer-as-library**:
Components ship library code (`src/`) and small static metadata; workflows hold the scripts that import the libraries and generate every workflow artifact, and own the data/weights they consume (ADR 0008). A component does not write into a workflow's folder at runtime; workflows regenerate artifacts by running their scripts.
_Avoid_: dependency, integration
