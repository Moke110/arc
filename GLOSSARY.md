# ARC

AI-driven drug repositioning for cardiac diseases. Predict the transcriptional perturbation a drug causes on cardiomyocytes, then rank drugs by predicted MIRI efficacy. The repository holds five external-source data/model categories (one directory per download source; see *External data provenance*), the Stage 1 drug-repositioning workflows (`stage1/`), Stage 2 model setup/selection (`stage2/`), and dev-tool files.

## Language

### Roadmap

**Stage**:
A program-level milestone of the AICM roadmap, with its own entry and exit conditions, spanning one or more workflows. A Stage is not a workflow *step* (one of the 11 numbered steps inside a workflow) and not a *component* (a reusable static library).
_Avoid_: phase, milestone, step

**Stage 1**:
The completed AICM stage that established XPert+DToxS MIRI drug repositioning: the four MIRI screens v1 (frozen baseline), v1.1 (DToxS-LINCS testbed), v2 (primary+iPSC), and v3 (all-CM, current frontier). "DToxS" names their shared primary-arm test-set anchor, not their only corpus.

**Stage 2**:
The AICM stage that selects the open-source whole-transcriptome drug-perturbation model that should feed future MIRI perturbome prediction. Setup (`stage2/01_setup/`) proved nine deployed non-preprint models load or train in their own environments; `stage2/02_dataset_corpus/` owns the five source-dataset corpora and their crosswalk tables (ADR 0033, amending ADR 0026); the cold-drug benchmark that compares the seven capable models is `stage2/03_models_cmbench/` (frozen pending a redo on the merged corpus; ADR 0027); `stage2/05_prnet_pretrain/` owns the CM-schema Tahoe view and the pretrain run on it (executed 2026-11-01; ADR 0031); `stage2/06_prnet_transfer/` is the concluded transfer study (ADR 0032). `04` is a deliberate gap: the retired Tahoe data step was folded into `02_dataset_corpus/`. Stage 2 evaluation excludes 978-landmark artifacts; PRnet's 978 demo is setup-only.
_Avoid_: phase; Stage 2 v1/v2 (the old versioned benchmark plan is archived)

**Stage 2 tree**:
The top-level `stage2/` directory that runs model setup and, later, model selection. It is distinct from a Stage 1 *repositioning workflow*; it produces no MIRI drug ranking in the setup stage. Preprint models are excluded from its setup set.

**Stage 3**:
The AICM stage that predicts the direction of change of one cardiomyocyte phenotype under one specified experimental condition from a drug-induced perturbome. Its first effort establishes the drug–phenotype data foundation; the perturbome will come from Stage 2 model predictions.

### Stage 3 phenotype data

**Phenotype**:
The biological indicator represented by a Stage 3 corpus target. MTT and CCK8 measurements both belong to `cell_viability`; their assay methods and source measurement scales remain identifiable.

**Condition**:
The detailed biological and experimental context of a Stage 3 phenotype measurement, including cell type and drug, hypoxia, and reoxygenation timings. Drug time denotes pre-HR treatment duration in an HR experiment or direct exposure duration without HR; drug persistence during HR remains separately described.

**Readout**:
The concrete numerical measurement provided by an experimental source, together with its unit, assay method, and raw or author-normalized representation. A readout is distinct from both the phenotype being measured and its endpoint label.

**Endpoint**:
The numeric corpus response class: +1 or −1 requires both the corresponding change-magnitude criterion and statistical significance; 0 is the rule's non-change class when the joint criterion fails on evaluable evidence. Class 0 does not itself assert demonstrated biological equivalence, and final corpus endpoints are numeric rather than an unknown class.

**Phenotype corpus sample**:
A source-dataset–drug–dose–endpoint unit within a particular condition and phenotype. Records from the same dataset with the same remaining key are grouped together, retaining all readout values rather than averaging them.

**Phenotype corpus**:
The quality-controlled collection of phenotype corpus samples, partitioned by fine condition and phenotype, with numeric endpoints, retained measurement sequences and provenance. It precedes the selection and aggregation performed by model data views.

**Readout sequence**:
The numerical source measurements retained by one merged phenotype corpus sample, preferring raw values and retaining author-normalized values separately. A normalized-only source retains its supplied values with that representation identified; per-value assay, unit, experiment and source references remain recoverable.

**Phenotype fold change**:
The ratio of a treated-group mean to the mean of its matched control on a compatible measurement scale, requiring a positive, valid denominator. It describes an experimental phenotype contrast rather than a transcriptomic fold change.

**BH family**:
The predeclared set of eligible drug–dose comparison tests sharing source dataset, full condition, phenotype, assay method and original experimental batch, whose p-values are adjusted together. Its membership is determined before selecting tests by effect magnitude or significance.

**Model data view**:
A selection of corpus conditions and samples, with explicitly chosen merging and aggregation rules, for a particular modeling task. These choices do not overwrite the underlying corpus measurements.

**Phenotype data**:
Experimental measurements linking pharmacological perturbations to cell phenotype readouts and their controls, conditions, and provenance. These are distinct from transcriptomic perturbation measurements and from model-predicted perturbomes.

**Phenotype observation**:
A drug-associated measurement of one cardiomyocyte readout under a documented experimental condition, with its experimental source and reference control retained. Multiple doses, times, cell models, and repeats must remain distinguishable.
_Avoid_: independent sample (without a counting definition)

**Qualified phenotype observation**:
A drug–readout observation that meets the declared experimental-condition, identity, and measurement-evidence criteria of a Stage 3 counting stratum. Qualification for counting does not by itself establish a reliable direction label for training.
_Avoid_: training sample (when only counting eligibility is established)

**Stage 3 drug identity**:
The evidence-resolved active parent of a single pharmacological intervention, retaining its known stereochemistry. The administered salt, hydrate, or other chemical form remains separately described; an alias or a different salt does not by itself create another drug.

**Phenotype counting stratum**:
The comparable target used for evidence counting in the wayfinder design, fixing a measurement definition, biological model, protocol, vehicle profile, and reference meaning while allowing declared dose/time axes. It is distinct from the finer condition partitions of the formal Stage 3 corpus.

**Phenotype combination**:
The wayfinder evidence-counting unit: one drug under one exact condition measured for one indicator within a phenotype counting stratum. It is distinct from a corpus sample whose identity also retains endpoint.

**Readout direction**:
The increase, non-change, or decrease classification relative to a specified control under the declared endpoint rule. Whether that direction implies benefit depends on the phenotype and condition.
_Avoid_: efficacy label (when only a measurement direction is known)

**Label-ready phenotype observation**:
A phenotype observation admitted as a numerically labeled corpus record under the frozen magnitude-and-significance rule, with its input and processing evidence retained. Point estimates or source-table coverage alone do not establish this status.

**Unresolved direction**:
A source-level observation for which the endpoint rule cannot yet be evaluated from the available evidence. It is not an allowed final corpus endpoint or an additional phenotype class.

### External data provenance

The repository has five external-source categories — `pert_data/`, `pert_model/`,
`disease_data/`, `drug_emb_model/`, and `pheno_data/` — holding 32 active 项 in total
(8 + 14 + 3 + 2 + 5), one per download source. ADR 0038 adds the independent
phenotype category to ADR 0030's historical four-category, 27-item inventory. Each 项 is verified offline from its tracked manifest by
`script/verify_all.sh`; its machine-readable record is the author log.

**External item (项)**:
One directory inside an external-source category that corresponds to exactly one download source — the unit at which provenance is recorded, named after the upstream product. Items are disjoint and each is self-contained.
_Avoid_: dataset (when meaning the unit), download

**Sub-item (项内成员)**:
One of the fixed members that make up an external item: its download scripts, its source tree, its declared local additions, its author log, and its manifest. The membership is fixed by the layout contract, not a matter of taste or convenience.
_Avoid_: subdirectory, component

**Source tree**:
The part of an external item that is byte-identical to what upstream shipped, kept in upstream's own layout and never edited in place — it is replaced only by re-downloading. Being untouched is what lets it be diffed against a fresh download.
_Avoid_: raw data, cache

**Local addition**:
Content inside an external item that did not come from upstream — generated here or extended here. Every local addition is declared individually with its origin and kept apart from the source tree, so the latter stays comparable to upstream.
_Avoid_: patch, modified download

**Author log**:
The machine-readable provenance record of one external item: what the item is, where each part came from, which version or commit, when it was retrieved, and whether it still verifies. It is the durable answer to "what is this, and is it intact?".
_Avoid_: metadata (unqualified), provenance notes

**Manifest**:
The per-file hash list over one external item's source tree, with which the item can be verified offline — no network and no original archive required — by recomputing the hashes.
_Avoid_: checksum list (when meaning the whole record), lockfile

**Execution area**:
The part of a Stage 2 step that holds the harness and the intermediates — code, configurations, run cards, tests, logs, checkpoints. Agents own it; it is not a result and no human reads it as a conclusion.
_Avoid_: code folder, scripts

**Output area**:
The part of a Stage 2 step that holds artifacts a later step reads or a human must know exists — frozen inputs, predictions, expression matrices. The conclusions a human opens (figures and result tables) sit above it, at step level.
_Avoid_: results (when meaning the whole area), artifacts

**Substep**:
An independently concluded part of a Stage 2 step, placed in its own `NN.N_<outcome>/` directory that carries the full step layout, including its own human-facing figure and table. The parent step then keeps only its own two documents and the substep directories.
_Avoid_: phase, sub-stage, chapter

### Categories

**Dev-tool files**:
Files that support the workspace itself rather than the science: `.gitignore`, `skills-lock.json`, `.agents/`, `docs/`, `reference/`, `presentation/`, `AGENTS.md`, `environment.yml`.
_Avoid_: tooling, meta

**Perturbation component**:
One of the three building blocks of prediction: a perturbation dataset (`pert_data/`), a drug embedding model (`drug_emb_model/`), or a perturbation prediction model (`pert_model/`). Components are pure static libraries (ADR 0008): `src/` code + small metadata only; downloaded data/weights live with the external item that provenance-tracks them (an item's `source/` tree, ADR 0030) or in the consuming workflow.
_Avoid_: module

**Repositioning workflow**:
A complete end-to-end drug-screening pipeline for one disease, bundling its scripts, prompts, data, and notebook. MIRI workflows live under `stage1/`: `stage1/xpert-fda-miri-v1/` is the canonical finalized MIRI screen (11 steps, frozen DToxS finetune data, untouched); `stage1/xpert-fda-miri-v1.1/` is its DToxS-LINCS testbed sibling (finetune data rebuilt from the LINCS LDS-1593 wells, per-plate normalized since ADR 0016, 11 steps — the agentic report step was removed); `stage1/xpert-fda-miri-v2/` is its DToxS-primary + DToxS-iPSC-CM finetune sibling (ADR 0015), aligned to v1.1's benchmark test sets with a single-cell encoding, per-plate normalized since ADR 0016, 11 steps; `stage1/xpert-fda-miri-v3/` is its all-CM finetune sibling (ADR 0017) training on all four usable CM datasets (primary + iPSC + gse244740 + gse262419) with train-only benchmark rows for the added datasets, 11 steps; `stage1/xpert-fda-miri-v0/` is the frozen legacy snapshot.
_Avoid_: pipeline (when meaning a workflow), project

**Disease data**:
Shared disease-specific inputs, reusable across future workflows: raw GEO acquisitions and per-GSE DEG lists under `disease_data/`, one external item per GEO series — upstream-derived files in the item's source tree, locally derived DEG lists in its declared local additions. Distinct from workflow artifacts because more than one workflow can consume it.
_Avoid_: input data, shared data

### Components

**Drug embedding**:
A learned vector representation of a drug: HG biological embedding (256-dim) and UniMol 3D features (122×514). Produced by embedding-model inference; consumed by finetune and predict. The finalized MIRI screen uses only the author-reference HG table (byte-identical) plus the precomputed UniMol array.
_Avoid_: drug feature

**Perturbome (扰动组)**:
The transcriptomic change Δ a drug induces in a cell relative to its matched control — scale-agnostic, the umbrella notion. Two scales in AICM: the **978 L1000-landmark** instance (Stage 1, see `Perturbed transcriptome`) and the **whole-transcriptome** instance (Stage 2 model output). The thesis spine is 模型→药效→闭环 over the perturbome.
_Avoid_: perturbation profile, drug signature

**Perturbed transcriptome**:
The **978-landmark instance** of the perturbome: the predicted expression of 978 L1000 landmark genes for a drug-treated cell. Produced by XPert inference; z-scored and scored in the workflow. Distinct from Stage 2's whole-transcriptome instance.

**Drug list**:
The set of drugs a workflow screens — for `xpert-fda-miri-v1`/`v1.1`/`v2`/`v3`, the 1,004 L1000-contained FDA drugs (exact canonical-SMILES match, `drug_idx < 8981`). Identity of each drug is its SMILES; the canonical mapping is `all_drugs_idx2smi.npy`.
_Avoid_: drug set

**Baseline transcriptome**:
The consensus control expression (978 genes) of the cell type under study, derived from the DToxS-primary control wells. Workflow-owned; the reference condition the perturbed transcriptome is measured against.

**Union finetune data**:
The workflow-built finetune h5ad of the CM datasets (step `03_cm_finetune_data`), used by the rebuilt workflows' finetune/benchmark/predict. HG-resolved compounds only, 978 landmarks, per-sample matched controls, precomputed dose/time indices. Since ADR 0016 every sample is **per-plate normalized**: its matched control is the mean of the control wells on its own `(dataset, plate)` and the 978-landmark matrix is quantile-aligned to the L1000 pretrain per `(dataset, plate)` group. v2's union = `dtoxs_primary.h5ad` + `dtoxs_ipsc.h5ad` (958 samples / 45 unique drugs); v3's union (ADR 0017) adds `gse244740.h5ad` (~2,424) + `gse262419.h5ad` (~936) for ~4,318 samples — gse244740-dominated (~56%), not gse262419 (its TempO-Seq panel is mostly environmental/toxicology compounds absent from the author table). The DToxS-primary arm reads the LINCS LDS-1593 well-level package and both arms carry real-concentration dose bins. The benchmark test sets are byte-identical to v1.1's fold test sets (v1.1 fold columns imported for the primary arm). Distinct from v1's frozen DToxS-only `stage1/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad` (dataset-wide normalized, the frozen historical anchor) and from v1.1's DToxS-LINCS finetune data.

**Per-plate normalization**:
The step-03 method (ADR 0016) that groups both the quantile alignment and the matched-control computation by `(dataset, plate)`: each sample's control baseline is the mean of the control wells on the same plate, and each plate's 978-landmark distribution is aligned to the L1000 pretrain independently. Applies to every dataset in the rebuilt finetune builds (v1.1/v2/v3); v1's frozen `stage1/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad` stays dataset-wide and is the documented historical anchor. Step 02's consensus baseline is *not* per-plate — it is a single dataset-wide scoring reference.
_Avoid_: per-batch normalization, plate-corrected

**DToxS-LINCS finetune data**:
The workflow-rebuilt DToxS finetune h5ad (662 samples / 45 drugs) that v1.1 trains on, built by its step `03_cm_finetune_data` from the LINCS LDS-1593 well-level package under `pert_data/LINCS-LDS-1593/source/Data/Level3` (released flat 2026-08-18) — a legitimate superset of v1's frozen 521/36 data (the raw GEO deposit lacks the non-KI perturbations). Schema-parity-gated against v1's frozen `stage1/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad` with an auditable delta report.

**LINCS LDS-1593**:
The DToxS center's well-level Level-3 deposit (LINCS LDG-1452), released flat at `pert_data/LINCS-LDS-1593/source/` on 2026-08-18 (the former `lincs-lds-1593/` wrapper removed; the GSE146095/GSE146096 raw GEO deposit was deleted the same day as a strict subset). Carries the full 956-well campaign (23 KIs + non-KI drugs + anthracyclines + antibody/protein wells) and the LINCS metadata suite (SMILES, batch codes). The v1-era package is frozen byte-identically at `pert_data/LINCS-LDS-1593/local/v1-frozen/` (incl. v1's frozen `dtoxs.h5ad`).

**Author drug embedding table**:
The author-provided L1000 HG embedding table (rows 0–8980, `drug_idx < 8981`), workflow-owned at `stage1/xpert-fda-miri-v1/input/HG_drug_embeddings.npy`. A compound is an **in-HG drug** when it resolves to a row of this table (canonical-SMILES or author `cmap_name` match); compounds resolving only to the FDA-extension rows (≥ 8981) or not at all are excluded from finetune and recorded in the step-03 manifest.
_Avoid_: HG-resolved, drug-set member

**Collapsed CM cell identity**:
The convention that every CM finetune sample carries one shared cell identity — `cell_iname="cardiomyocyte"`, `cell_idx=0`, `tissue_idx=0` (`num_cell_id 1`) — so the model never distinguishes donors or cell lines. Biological (donor, drug) identity is not model-visible; the benchmark fold gates assert well-identity (no identical well in train and test) and drug identity (cold-drug test drugs excluded from every finetune dataset's train). A model-schema convention, not a statement about data modality: every CM source is a **Bulk profile**.
_Avoid_: single-cell CM encoding (reads as data modality), cardio cell, per-donor encoding

**Finetune concentration**:
The per-drug real µM concentration used to bin each finetune sample's dose via the author `assign_pert_dose_corrected` scheme. Sources: the DToxS-primary paper Table S1 (23 KIs), the iPSC-CM series matrix (shared non-KIs, verified 22/23 match), and literature-curated clinical Cmax with citations for the remaining 8 primary-only drugs (SB-431542 uses its standard in-vitro assay concentration — it has no clinical Cmax). This amends the "no external database" provenance rule for the dose field only: concentrations must be cited.
_Avoid_: single-dose convention, dose 0

### Transcriptome scale

**Bulk profile**:
The transcriptome of one physical sequencing library pooled from many cells — one library per culture well. All four CM sources are this kind (two whole-transcriptome RNA-seq, one ScreenSeq 3′ UMI plate RNA-seq, one TempO-Seq targeted panel).
_Avoid_: bulk data, bulk expression, bulk (unqualified)

**Pseudobulk profile**:
Counts summed across the cells of one cell line inside one well (or one plate-matched group) of a single-cell atlas — an analysis-level aggregate, not a physical library. The bulk form of Tahoe-100M is this kind, never a Bulk profile.
_Avoid_: bulk profile, pseudo bulk, pooled single cell

**Well**:
One culture well holding one perturbation condition: the row key of the CM corpora and of the Tahoe pseudobulk corpus. Tahoe-100M calls it a `sample` (96 wells per plate, 1,344 wells in 14 plates), so the two vocabularies must be translated explicitly.
_Avoid_: sample (when meaning well), culture

### Workflow step layouts

All four MIRI screens (v1, v1.1, v2, v3) share the same 11-step layout: `01_drug_embed`, `02_cm_baseline`, `03_cm_finetune_data`, `04_xpert_pretrain`, `05_xpert_finetune_bench`, `06_xpert_finetune`, `07_xpert_predict`, `08_miri_phenotype`, `09_miri_efficacy_scores`, `10_filter`, `11_integrate`. They differ only in the finetune-data composition (v1's frozen DToxS / v1.1's DToxS-LINCS / v2's primary+iPSC union / v3's all-four-CM union). The agentic report step is **no longer a numbered step** in any workflow (removed per spec 1); v1 keeps its archived, non-numbered `report/` stage.

**Non-numbered stage**:
A workflow subdirectory that is neither one of the numbered screen steps nor a *component*: leading-underscore material (`_lib/`, v2's `_other/`) or a workflow's archived stage (v1's `report/`). Signals "not part of the screen pipeline".
_Avoid_: extra step, utility folder

**Limitations analysis**:
A read-only diagnostic step that quantifies a screen's structural bounds — landmark coverage, critical-pathway overlap, screened-drug embedding structure — without changing any screen artifact. v2's is `13_limits/`.
_Avoid_: weakness report, ablation study

**Pretrain unit**:
The workflow-owned author pretrained XPert weight + pretrain YAML (`04_xpert_pretrain/` in all of v1/v1.1/v2/v3); the L1000 pretrain h5ads live in `pert_data/XPert-L1000/source/data/` (item source tree). Both pretrain and finetune YAML recipes are owned by the workflow.

**MIRI efficacy score**:
A score of how strongly a drug's predicted perturbation reverses the MIRI transcriptional disease signature. Computed by three methods (WCS+RGES, Zhang/ssCMap, pathway RGES) yielding four z-components.
_Avoid_: drug score

**Cardiotoxicity score**:
Deferred from v1. A score of a drug's risk of cardiac harm (e.g., hERG/QT).
_Avoid_: safety score

**MIRI reference DEG list**:
The integrated, cross-dataset MIRI differential-expression reference (`MIRI_up-down.json`, `MIRI_rank.json`) built from the per-GSE DEG lists. Distinct from per-GSE DEG lists and from the curated `MIRI_pathways.json`.

**Per-GSE DEG list**:
The DEG list derived from a single GEO study, kept alongside its acquisition data in `disease_data/<GSE>/local/deg/`. The raw material for the MIRI reference DEG list.

**Final drug report**:
The agentic bilingual (EN/CN) per-drug fact-sheet report, driven by `prompt.md` and assembled from resumable per-drug agent JSONs. Exists only in v1 as an **archived, non-numbered `report/` stage**; the rebuilt workflows (v1.1/v2/v3) have **no report stage** — they end at integrate.
_Avoid_: drug research report

### Stage 2 setup

**Setup stage**:
The current Stage 2 stage: each deployed non-preprint model is loaded or trained on its native asset and produces one prediction artifact, a run log, and a run card. It answers "can this model run here?" and deliberately does not compare models or compute a benchmark metric.
_Avoid_: Phase 0 (plan shorthand only), smoke test (when meaning the whole stage)

**Run card**:
The one-page record for a model's setup: source commit, weights, data, pinned environment, command, expected artifact, and caveats. Lives at `stage2/01_setup/execution/<model>/README.md`.
_Avoid_: model README, config

**Native asset**:
The checkpoint and dataset a model's authors ship or document for that model. A native-asset run may use a non-chemical dataset (CPA/Kang) or a landmark subset (PRnet/LINCS); that limitation is recorded on the run card.

**Setup-only artifact**:
An output permitted only to prove that a model runs, not to support a Stage 2 claim. PRnet's 978-gene LINCS demo is the canonical example; it is quarantined under `stage2/01_setup/execution/prnet/` and never read by Stage 2 evaluation.

**Selection stage**:
The deferred Stage 2 stage after setup: choose which model(s) to adapt to CM data and define the evaluation. Its metric, folds, and gene space are not fixed by the setup plan.

### Stage 2 dataset corpus

**Dataset corpus**:
The information-maximal h5ad substrate between a raw downloaded dataset and a model's training input: one file per source dataset, raw counts in `X`, the source's **native** gene axis, every raw well including vehicle controls and technical blanks, and no model-facing filter, no integer encoding and no Δ. Lives in `stage2/02_dataset_corpus/output/corpus/<version>/`; the five corpora are `dtoxs_primary`, `dtoxs_ipsc`, `sctox`, `tenv` and `tahoe`.
_Avoid_: CM perturbation corpus (retired name for the matrix-free gene/drug universe), union finetune data, finetune corpus, training set, dataset bundle

**Corpus coverage**:
The README record of what a corpus contains, as five counts: distinct cell lines, distinct chemicals, distinct doses, distinct times, and distinct `(cell line, chemical, dose)` combinations at that corpus's own time point; vehicle controls and technical blanks are counted separately and never as chemicals.
_Avoid_: corpus stats, dataset summary

**CM gene axis**:
The frozen 18,164-symbol four-source intersection named by latest HGNC symbol — the shared axis a CM benchmark or CM-model view projects onto, not a corpus gene universe (each corpus keeps its source's native axis).
_Avoid_: corpus gene set, CM universe

**Dataset id**:
The short Stage 2 name of a source dataset — `dtoxs_primary`, `dtoxs_ipsc`, `sctox`, `tenv`, `tahoe` — carried in each h5ad's `dataset_id` column. The GEO accession or atlas name stays the provenance anchor on the raw folder; historical records that predate the 2026-11 rename keep `screen_ctox`, `temp_env` and `dtoxs_cm`.
_Avoid_: source name, dataset name, accession

**View**:
A derived, model-facing h5ad that a consuming step builds from one or more corpora (one pharmacology tier, or a projection onto the CM gene axis); it is reproducible and disposable, and it is owned by the consumer that defines its rule — the corpus step never ships one.
_Avoid_: subset, variant, dataset bundle

**Pretrain view**:
A view that projects one corpus's gene axis and column conventions onto another corpus's (today: `tahoe` onto the CM gene axis), so a model pretrained on the source can be adapted to the target without a gene-space remap. Produced by the pretraining step that consumes it (`stage2/05_prnet_pretrain/`), not by the corpus step.
_Avoid_: pretrain corpus (for a view — the corpus is the source, never the projection), pretrain dataset, projection (unqualified)

**Dataset bundle** (retired 2026-11-04):
The former numbered view build under `stage2/02_cm_data/output/datasets-<n>/`; replaced by the dataset corpus plus consumer-owned views (ADR 0033).
_Avoid_: dataset version, output, export

**Frozen CM evidence (v1)**:
The archived five 18,164-gene CM views of the retired bundle, kept read-only because a frozen benchmark and a concluded transfer study name them as input evidence; nothing new is built from it.
_Avoid_: legacy bundle, old views

**Pharmacology tier**:
The named compound-inclusion rule of a view: `atc-moa` (ATC or ChEMBL mechanism), `atc-moa-mesh` (plus PubChem MeSH classes), or `all-smiles` (any resolvable small molecule).
_Avoid_: definition 1/3, filter level, drug class

**Drug key**:
The desalted parent InChIKey (stereochemistry retained) that identifies a drug across names, salt forms, and datasets.
_Avoid_: drug id, compound name, drug_idx

**Matched control group**:
The `dataset::plate` group of vehicle wells a perturbed sample is compared against.
_Avoid_: baseline, control set, batch control

**Manual pharmacology supplement**:
A cited, per-compound completion of the pharmacology annotation, allowed only on the curated drug panels (`dtoxs`, `gse244740`) and never on the `gse262419` tiers.
_Avoid_: manual override, curation fix

**Pharmacology-qualified compound (药理资格合格化合物)**:
A human-use approved or clinical-stage drug, or a research tool compound with a documented pharmacological mechanism, admitted by evidence review for the formal PRnet CM study. This qualification is separate from an automatic pharmacology tier and from experimental quality or transcriptional response strength.
_Avoid_: ATC-positive drug, trainable drug (when only pharmacology has been reviewed)

**CM consensus control (CM 共识对照)**:
A fixed reference expression profile derived from eligible training CM vehicle controls, representing the CM background on which novel-drug effects are predicted. It is distinct from a treated observation's experimentally matched control and requires no new control measurement at prediction time.
_Avoid_: matched control, Stage 1 baseline transcriptome


### Stage 2 benchmark

**Tenv benchmark sample**:
One model's training-and-test result under one fixed fold and one fixed training
seed in the redesigned tenv comparison. Three folds crossed with three seeds
yield nine run samples per model; these are not nine independent biological
samples.
_Avoid_: well, drug–dose group, independent replicate

**Cold-drug fold**:
One of five Stage 2 benchmark splits in which a drug is never observed in the fold's train or validation set and is scored exactly once as test; every drug is test exactly once and valid exactly once across the five folds (6:2:2 = 276/92/92 drugs).
_Avoid_: cold-start split, holdout (unqualified)

**CM-unseen drug (CM 未见药物)**:
A drug whose cardiomyocyte perturbation observations are withheld from model fitting across all CM sources used for adaptation. The drug may have been observed in non-CM pretraining, including Tahoe; this is the prediction setting for drugs without available CM wet-lab responses.
_Avoid_: globally unseen drug, unseen drug (without specifying the CM scope)

**Δ space**:
The benchmark's target: `log1p_cp10k(treated) − mean(log1p_cp10k(same-plate DMSO))`, in the 18,164-gene corpus order. The model emits per-well Δ; the runner averages within a group.
_Avoid_: fold change, expression delta (unqualified)

**Group**:
One `(drug_key, dose)` pair in the benchmark — 1,380 groups, each the mean of its 3–6 wells. AUPRC is macro-averaged over groups.
_Avoid_: condition (unqualified), sample

**Cold-drug-capable model**:
A model whose drug input is derived from chemical structure (fingerprint, RDKit descriptor, or a structure encoder), so it can score a drug it never saw in training. Stage 2's admission gate.
_Avoid_: warm-drug model, structure model

**Entrant set**:
The seven cold-drug-capable models in the Stage 2 benchmark: PRnet, CPA, chemCPA, biolord, PerturbNet, Squidiff, CMonge. scVIDR and DECANT are documented exclusions, not rows.
_Avoid_: candidate list, all models

**Ideal-fit label**:
The `ideal_fit` report column: `full` (author architecture is used as intended under the contract), `partial` (a declared omission or transductive control, e.g. chemCPA's skipped 978 pretrain, biolord's transductive prediction, PerturbNet's missing control transcriptome), or `none` (a CPU baseline). One leaderboard; the label is not a second table.
_Avoid_: ideal model, tier

**Benchmark run**:
One `(model, fold)` unit in `stage2/03_models_cmbench/03.1_cold_drug_18k/output/units/`: fit on the training drugs, select the checkpoint on the cold validation drugs, predict the test wells, and write the per-unit artifacts and `status.json`.
_Avoid_: experiment, trial


### Stage 1 evaluation

**Top-k true set**:
The ground-truth label convention of the Stage 1 XPert benchmark — the step-05 official folds and the `warmstart_pr/` PR/AUPRC analysis: per sample, the true set is the top-k up genes and the bottom-k down genes of that sample's true Δ (k = 100 in the 978-landmark space), every other gene being a negative in both directions. A fixed-truncation (ranking) definition with no significance test, so each direction's prevalence is exactly k/978 and that prevalence is the random AUPRC baseline.
_Avoid_: DEG set (unqualified), differential set

Historical Stage 2 v1 benchmarking terms (direction-aware AUPRC, unified test set, data-ablation benchmark, dtoxs-primary benchmark) are archived with `stage2/_archive/cm_bench_v1/` and are not active vocabulary.

### Experimental validation

**Wet-lab validation (Stage 1)**:
The first experimental check of Stage-1 drug candidates, run outside the repository and recorded at `docs/wetlab-validation-stage1.md` (raw data in `wet_data/`, gitignored): cardiomyocyte CCK8 viability (safety; hypoxia/reoxygenation efficacy) and a mouse IR model (%EF/%FS/infarct, with a non-IR safety arm). The two candidates are trimethoprim (MP) and mirtazapine (MT).
_Avoid_: in-vivo validation, experiment (unqualified)

**MP / MT**:
Two-letter shorthand for the Stage-1 wet-lab-validation candidates — **MP** = trimethoprim (甲氧苄啶), **MT** = mirtazapine (米氮平). Used only in the wet-lab-validation record and discussions of it; the screens themselves use full drug names.
_Avoid_: reusing MP/MT for other compounds


### Boundary

**Producer-as-library**:
Components ship library code (`src/`) and small static metadata; workflows hold the scripts that import the libraries and generate every workflow artifact, and own the data/weights they consume (ADR 0008). A component does not write into a workflow's folder at runtime; workflows regenerate artifacts by running their scripts.
_Avoid_: dependency, integration
