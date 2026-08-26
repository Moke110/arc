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
A complete end-to-end drug-screening pipeline for one disease, bundling its scripts, prompts, data, and notebook. MIRI workflows live under `workflow/`: `workflow/xpert-fda-miri-v1/` is the canonical finalized MIRI screen (11 steps, frozen DToxS finetune data); `workflow/xpert-fda-miri-v1.1/` is its DToxS-LINCS testbed sibling (12 steps, finetune data rebuilt from the LINCS LDS-1593 wells, v1 untouched); `workflow/xpert-fda-miri-v2/` is its DToxS-primary + DToxS-iPSC-CM finetune sibling (ADR 0015), aligned to v1.1's benchmark test sets with a single-cell encoding; `workflow/xpert-fda-miri-v0/` is the frozen legacy snapshot.
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
The set of drugs a workflow screens — for `xpert-fda-miri-v1`/`v2`, the 1,004 L1000-contained FDA drugs (exact canonical-SMILES match, `drug_idx < 8981`). Identity of each drug is its SMILES; the canonical mapping is `all_drugs_idx2smi.npy`.
_Avoid_: drug set

**Baseline transcriptome**:
The consensus control expression (978 genes) of the cell type under study, derived from the DToxS-primary control wells. Workflow-owned; the reference condition the perturbed transcriptome is measured against.

**Union finetune data**:
The workflow-built finetune h5ad of the CM datasets (step `03_cm_finetune_data`), used by v2's finetune/benchmark/predict. HG-resolved compounds only, 978 landmarks, per-sample matched controls, precomputed dose/time indices. Built by concatenating the per-dataset h5ads (`dtoxs_primary.h5ad` + `dtoxs_ipsc.h5ad` = 958 samples / 46 unique drugs); the DToxS-primary arm reads the LINCS LDS-1593 well-level package and both arms carry real-concentration dose bins. The benchmark test sets are byte-identical to v1.1's fold test sets (v1.1 fold columns imported for the primary arm). Distinct from v1's frozen DToxS-only `input/dtoxs.h5ad` and from v1.1's DToxS-LINCS finetune data.

**DToxS-LINCS finetune data**:
The workflow-rebuilt DToxS finetune h5ad (662 samples / 45 drugs) that v1.1 trains on, built by its step `03_cm_finetune_data` from the LINCS LDS-1593 well-level package under `pert_data/dtoxs_primary/Data/Level3` (released flat 2026-08-18) — a legitimate superset of v1's frozen 521/36 data (the raw GEO deposit lacks the non-KI perturbations). Schema-parity-gated against v1's frozen `input/dtoxs.h5ad` with an auditable delta report.

**LINCS LDS-1593**:
The DToxS center's well-level Level-3 deposit (LINCS LDG-1452), released flat at `pert_data/dtoxs_primary/` on 2026-08-18 (the former `lincs-lds-1593/` wrapper removed; the GSE146095/GSE146096 raw GEO deposit was deleted the same day as a strict subset). Carries the full 956-well campaign (23 KIs + non-KI drugs + anthracyclines + antibody/protein wells) and the LINCS metadata suite (SMILES, batch codes). The v1-era package is frozen byte-identically at `pert_data/dtoxs_primary/v1-frozen/` (incl. v1's frozen `dtoxs.h5ad`).

**Author drug embedding table**:
The author-provided L1000 HG embedding table (rows 0–8980, `drug_idx < 8981`), workflow-owned at `workflow/xpert-fda-miri-v1/input/HG_drug_embeddings.npy`. A compound is an **in-HG drug** when it resolves to a row of this table (canonical-SMILES or author `cmap_name` match); compounds resolving only to the FDA-extension rows (≥ 8981) or not at all are excluded from finetune and recorded in the step-03 manifest.
_Avoid_: HG-resolved, drug-set member

**Single-cell CM encoding**:
The convention that every CM finetune sample carries one shared cell identity — `cell_iname="cardiomyocyte"`, `cell_idx=0`, `tissue_idx=0` (`num_cell_id 1`) — so the model never distinguishes donors or cell lines. Biological (donor, drug) identity is not model-visible; the benchmark fold gates assert well-identity (no identical well in train and test) and drug identity (cold-drug test drugs excluded from both datasets' train).
_Avoid_: cardio cell, per-donor encoding

**Finetune concentration**:
The per-drug real µM concentration used to bin each finetune sample's dose via the author `assign_pert_dose_corrected` scheme. Sources: the DToxS-primary paper Table S1 (23 KIs), the iPSC-CM series matrix (shared non-KIs, verified 22/23 match), and literature-curated clinical Cmax with citations for the remaining 8 primary-only drugs (SB-431542 uses its standard in-vitro assay concentration — it has no clinical Cmax). This amends the "no external database" provenance rule for the dose field only: concentrations must be cited.
_Avoid_: single-dose convention, dose 0

### Workflow stages (finalized step layouts)

The v1 MIRI screen reads top-down: `01_drug_embed`, `02_cm_baseline`, `03_xpert_pretrain`, `04_xpert_finetune_bench`, `05_xpert_finetune`, `06_xpert_predict`, `07_miri_phenotype`, `08_miri_efficacy_scores`, `09_filter`, `10_integrate`, `11_report`. The v1.1 and v2 siblings insert `03_cm_finetune_data` (finetune-data prep) and shift the rest to 12 steps: `01_drug_embed`, `02_cm_baseline`, `03_cm_finetune_data`, `04_xpert_pretrain`, `05_xpert_finetune_bench`, `06_xpert_finetune`, `07_xpert_predict`, `08_miri_phenotype`, `09_miri_efficacy_scores`, `10_filter`, `11_integrate`, `12_report`. v1.1 and v2 share this layout and differ only in the finetune-data composition (DToxS-LINCS only vs union).

**Pretrain unit**:
The workflow-owned author pretrained XPert weight + pretrain YAML (`04_xpert_pretrain/` in v2, `03_xpert_pretrain/` in v1); the L1000 pretrain h5ads live in `pert_data/l1000-xpert/` (dataset home). Both pretrain and finetune YAML recipes are owned by the workflow.

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
The step-11 (v1) / step-12 (v2) agentic bilingual (EN/CN) per-drug fact-sheet report over the `09_filter`/`10_filter` survivors with integrated-score percentile > 0.5, driven by `prompt.md` and assembled from resumable per-drug agent JSONs.
_Avoid_: drug research report

### Boundary

**Producer-as-library**:
Components ship library code (`src/`) and small static metadata; workflows hold the scripts that import the libraries and generate every workflow artifact, and own the data/weights they consume (ADR 0008). A component does not write into a workflow's folder at runtime; workflows regenerate artifacts by running their scripts.
_Avoid_: dependency, integration
