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
The AICM stage that selects an open-source whole-transcriptome drug-perturbation model for MIRI perturbome prediction. It separates model setup, the source-dataset corpus, and model comparison on held-out perturbations.

**Stage 2 tree**:
The top-level `stage2/` directory that runs model setup and, later, model selection. It is distinct from a Stage 1 *repositioning workflow*; it produces no MIRI drug ranking in the setup stage. Preprint models are excluded from its setup set.

**Stage 3**:
The AICM stage that predicts the direction of change of one cardiomyocyte phenotype under one specified experimental condition from a drug-induced perturbome. Its first effort establishes the drug–phenotype data foundation; the perturbome will come from Stage 2 model predictions.

### Stage 2 chemical representation

**ECFP4**:
The Stage 2 desalted-parent circular fingerprint with radius 2, atom connectivity
invariants, 2048 binary bits, and chirality enabled. RDKit's API is
`GetMorganGenerator`. The upstream PRnet fingerprint uses feature invariants
and 1024 bits, so it is a distinct representation. Sealed Stage 2 generations
record the project fingerprint as `morgan.npy` and `morgan`.

**ECFP4C**:
The name for the existing chirality-aware Stage 2 ECFP4 in the drug
representation comparison. The suffix C denotes chirality across the comparison’s
ECFP4C, ECFP6C, FCFP4C and FCFP6C arms.

**ECFP4_MLP**:
The Stage 2 fingerprint-only treated-expression baseline. The model takes the
2048 ECFP4 bits, applies a 128-unit ReLU hidden layer, and predicts 5000 genes
through a linear output. Sealed Stage 2 02 artifacts call it `MLP_Morgan`.

### Stage 3 phenotype data

**Stage 3 cold-drug evaluation**:
Evaluation on desalted drug parents absent from training, keeping all aliases and doses of each parent together.

**Stage 3 cold-dose evaluation**:
Evaluation at a dose level absent from training, for drugs observed at other training doses; interior held doses assess interpolation and boundary held doses assess extrapolation.

**Stage 3 warm-start evaluation**:
Evaluation on unseen drug-parent–dose combinations for which the parent and dose each occur in training. All source records for the same combination share one partition.

**Stage 3 Morgan baseline**:
A structure-based predictor of cardiomyocyte phenotype direction, serving as the comparator for a future predictor using drug-induced transcriptomic perturbation Δ. Both predictors target the same declared phenotype contrasts, with a separate predictor for each experimental-condition–phenotype pair.

**Phenotype**:
The biological indicator represented by a Stage 3 corpus target. MTT, CCK8, CellTiter-Glo (CTG) and imaging-based viability measurements belong to `cell_viability`; their assay methods and source measurement scales remain identifiable.

**Condition**:
The experimental state of a Stage 3 phenotype measurement, such as baseline or hypoxia–reoxygenation. Detailed context, including cell type and treatment, hypoxia and reoxygenation timings, further describes that state.

**Baseline**:
The source's normoxic reference condition for chemical perturbation, including its documented culture and exposure schedule. It includes treated measurements and does not mean vehicle-only controls.

**Hypoxia–reoxygenation (H-R)**:
The experimental condition in which a hypoxic interval is followed by restoration of oxygen. Its interpretation retains the source's hypoxia, reoxygenation and chemical-treatment schedule.

**Readout**:
The source-supplied measurement evidence, including individual values or reported group statistics such as mean, SD, n and p-value. Its interpretation includes the assay, units and raw or author-normalized scale.

**Endpoint**:
A derived response class describing a phenotype contrast under a declared labeling rule. An endpoint and its statistical test are not prerequisites for a Stage 3 corpus sample (ADR 0043).

**Phenotype corpus sample**:
A retained measurement record interpreted according to the source underlying its condition–phenotype corpus. Its representation follows the selected source; physically separate samples are not merged merely because their values match (ADR 0044).

**Phenotype corpus**:
The collection of source measurements for a chosen condition and phenotype, first built around its best available source and later extended with similar sources converted to its established representation. Source measurements retain experimental context and provenance without requiring significance tests or endpoint labels (ADR 0044).

**Replicate index**:
An index distinguishing physically separate samples whose other attributes apart from readout are identical. It does not imply biological independence.

**Source annotation**:
The description of one phenotype source’s experimental design, control types and readout meaning, with measurement-specific context when these vary within a source.

**Readout sequence**:
The collection of source measurement values retained together under the historical merged-sample design of ADRs 0039–0040. The current corpus preserves physically separate samples rather than merging them into such a sequence (ADR 0043).

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

**Stage 3 desalted parent**:
The chemical structure of a Stage 3 intervention after removing counterions and solvates and normalizing removable protonation, while retaining stereochemistry and permanent charges. Distinct source treatments may share a parent without becoming the same measurement.

**Phenotype counting stratum**:
The comparable target used for evidence counting in the wayfinder design, fixing a measurement definition, biological model, protocol, vehicle profile, and reference meaning while allowing declared dose/time axes. It is distinct from the finer condition partitions of the formal Stage 3 corpus.

**Phenotype combination**:
The wayfinder evidence-counting unit: one drug under one exact condition measured for one indicator within a phenotype counting stratum. It is distinct from the source-measurement identity of the current corpus sample.

**Readout direction**:
The increase, non-change, or decrease classification relative to a specified control under the declared endpoint rule. Whether that direction implies benefit depends on the phenotype and condition.
_Avoid_: efficacy label (when only a measurement direction is known)

**Magnitude-based direction**:
A decrease (-1), within-band (0), or increase (+1) label based on a declared effect-size band around a reference. It does not establish statistical significance or equivalence; a reference adopted as a working assumption remains explicit.

**Label-ready phenotype observation**:
A phenotype observation admitted as a numerically labeled corpus record under the frozen magnitude-and-significance rule, with its input and processing evidence retained. Point estimates or source-table coverage alone do not establish this status.

**Unresolved direction**:
A source-level observation for which the endpoint rule cannot yet be evaluated from the available evidence. It does not prevent retention of the source measurement in the phenotype corpus.

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
The build and audit machinery of a Stage 2 step, including its disposable working intermediates. It is distinct from reusable accepted products and the evidence establishing their acceptance.
_Avoid_: code folder, scripts

**Output area**:
The minimal set of accepted products a Stage 2 step supplies for reuse by other steps, including data, public evaluation interfaces and reusable trained models. Report figures and audit evidence are separate.
_Avoid_: results (when meaning the whole area), artifacts

**Stage 2 report figure**:
A human-facing visual conclusion whose complete plotted values are retained separately, allowing the same figure to be reproduced independently of model execution and scientific recalculation.

**Stage 2 acceptance evidence**:
The retained records establishing a product’s identity, provenance and successful scientific and downstream-use checks. Failure diagnoses and migration mappings preserve the interpretation of historical records.

**Substep**:
An independently concluded part of a Stage 2 step, placed in its own `NN.N_<outcome>/` directory that carries the full step layout, including its own report figures, plotted data and audit evidence. The parent step then keeps only its own two documents and the substep directories.
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

### Historical Stage 2 setup

**Setup stage**:
The retired Stage 2 setup stage: each admitted non-preprint model is loaded or trained on its native asset and produces one prediction artifact, a run log, and a run card. It answers "can this model run here?" and deliberately does not compare models or compute a benchmark metric.
_Avoid_: Phase 0 (plan shorthand only), smoke test (when meaning the whole stage)

**Run card**:
The one-page record for a model's setup: source commit, weights, data, pinned environment, command, expected artifact, and caveats. Historical run cards from the removed setup tree are retained in Git history.
_Avoid_: model README, config

**Native asset**:
The checkpoint and dataset a model's authors ship or document for that model. A native-asset run may use a non-chemical dataset (CPA/Kang) or a landmark subset (PRnet/LINCS); that limitation is recorded on the run card.

**Setup-only artifact**:
An output permitted only to prove that a model runs, not to support a Stage 2 claim. PRnet's 978-gene LINCS demo is the canonical example; its historical setup evidence is retained in Git history and is outside the current bases boundary.

**Selection stage**:
The historical proposal for a Stage 2 stage after setup: choose which model(s) to adapt to CM data and define the evaluation. Its metric, folds, and gene space are not fixed by the setup plan.

### Stage 2 dataset corpus

**Bases**:
The shipped Stage 2 foundation at `stage2/01_bases/`: frozen gene identity, raw counts, the integer corpus, the model-facing views and the released gene axes, plus disease and cardiotox DEGs, the per-disease fusion and the single public reader/metric API. It publishes `raw_counts`, `corpus`, `disease-deg`, `disease-fusion`, `gene-axis`, `gene-identity`, `palette`, `tenv-hvg`, `view-5k`, `api` and `fold_n_seed`. The retired measured truth, successor envelope and training corpora are not part of it; model setup and benchmark experiments are outside its scope.

**Raw-count preparation**:
Counts derived from immutable source expression before normalization or treated-well selection. Current bases retains tenv on its native gene axis and Tahoe projected onto the frozen CM gene axis, with controls retained.

**Corpus** (`output/corpus/`):
The integer raw-count substrate on the frozen CM gene axis: one row per treated or vehicle profile, with identity, dose, role and matched-control group, and with no normalization, no delta and no gene subsetting. Every gene-axis view and every model view is built from it; the former normalized-corpus meaning is retired.
_Avoid_: training corpus (the retired normalized collection), dataset corpus (the ADR 0033 native-axis per-source release)

**Dataset corpus**:
The information-maximal h5ad substrate between a raw downloaded dataset and a model's training input: one file per source dataset, raw counts in `X`, the source's **native** gene axis, every raw well including vehicle controls and technical blanks, and no model-facing filter, no integer encoding and no Δ. The historical five-source release comprised `dtoxs_primary`, `dtoxs_ipsc`, `sctox`, `tenv` and `tahoe`; its tracked evidence is retained in Git history.
_Avoid_: CM perturbation corpus (retired name for the matrix-free gene/drug universe), union finetune data, finetune corpus, training set, dataset bundle

**Corpus coverage**:
The README record of what a corpus contains, as five counts: distinct cell lines, distinct chemicals, distinct doses, distinct times, and distinct `(cell line, chemical, dose)` combinations at that corpus's own time point; vehicle controls and technical blanks are counted separately and never as chemicals.
_Avoid_: corpus stats, dataset summary

**CM gene axis**:
The frozen 18,164-symbol four-source intersection named by latest HGNC symbol — the shared axis a CM benchmark or CM-model view projects onto, not a corpus gene universe (each corpus keeps its source's native axis).
_Avoid_: corpus gene set, CM universe

**DEG (gene axis, v2)**:
A gene of one disease source dataset that passes `padj < 0.05` **and** `|log2FC| ≥ 2` (FC ≥ 4) in at least one of that dataset's contrasts; the multi-contrast collapse keeps the maximum-effect contrast. It replaces the Stage 2 01_bases v1 meaning, in which every measured gene carried a within-dataset rank. The tenv-HVG ranking is not a DEG source.
_Avoid_: significant gene (unqualified), hit, DE gene

**Source tier** (gene axis):
The evidence class of a disease source: **tier 1** human atria/ventricle with cardiomyocyte specificity (CM pseudobulk or CM-nuclear bulk); **tier 2** hiPSC/hESC-CM cells, organoids or engineered tissue (clinical nHF/ICM/AF admit only non-genetic models); **tier 3** human atria/ventricle bulk tissue without cardiomyocyte specificity. The ADR 0073 release collects all qualifying tiers without a cumulative DEG-count stop; ADR 0072's stopping rule is historical.
_Avoid_: priority (the retired 0066 name), level

**Highest DEG source tier** (`source_tier_highest`):
For one disease and one DEG, the smallest numeric tier among admitted informative sources in which that gene passes the DEG definition. A gene that is a DEG in tier 1 and tier 2 has highest tier 1; mere measurement is not DEG evidence.
_Avoid_: maximum tier number, best measured tier

**DEG1000** (ADR 0073 accepted design):
Per disease, the first 1,000 genes of its complete DEG union ordered by highest DEG source tier ascending, within-disease rank sum ascending, then frozen gene index ascending; a smaller union is retained whole. This replaces DEG2000 in the successor disease contribution to the gene axis.
_Avoid_: per-source top 1000, disease DEG signature

**f1000** (ADR 0073 accepted design):
The number of eligible diseases whose DEG1000 contains a gene; the disease-side replication count in the successor axis selection. It excludes tenv-HVG and diseases below the full-union eligibility gate.
_Avoid_: quality score, continuous rank, unqualified disease frequency

**DEG2000 and f2000** (historical release vocabulary):
The per-disease top-2,000 DEG set under the former disease-rank ordering and the count of eligible diseases containing a gene in that set. Retain these only when describing historical ADR 0071/0072 measured results; ADR 0073 uses DEG1000/f1000.
_Avoid_: DEG2000 or f2000 for the ADR 0073 successor

**Dataset id**:
The short Stage 2 name of a source dataset — `dtoxs_primary`, `dtoxs_ipsc`, `sctox`, `tenv`, `tahoe` — carried in each h5ad's `dataset_id` column. The GEO accession or atlas name stays the provenance anchor on the raw folder; historical records that predate the 2026-11 rename keep `screen_ctox`, `temp_env` and `dtoxs_cm`.
_Avoid_: source name, dataset name, accession

**View**:
The model-facing expression product derived from the integer corpus on a declared gene axis, carrying normalized expression, matched-control expression and their difference. ADR 0077 publishes only paired 5k expression with axis-local normalization; the ADR 0073 four-view full-axis release and existing downstream experiments remain sealed history.
_Avoid_: subset, variant, dataset bundle

**Pretrain view**:
A view that projects one corpus's gene axis and column conventions onto another corpus's (today: `tahoe` onto the CM gene axis), so a model pretrained on the source can be adapted to the target without a gene-space remap. The pretraining step that consumes the view owns it.
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

**Tenv highly variable genes (tenv HVG)**:
Genes on the frozen CM gene axis ranked by per-well expression dispersion across the **4,206 treated tenv wells** — the variability source of the gene axis, with no upregulation or downregulation direction. The train-only fold rankings in `output/tenv-hvg/` are a different population and are never an axis input.
_Avoid_: tenv HVG (unqualified, when the train-only model-selection ranking is meant)

**Stage 2 sample identity**:
The stable identity of one retained perturbation observation, distinguishing separate source wells and, for Tahoe pseudobulk, the cell line within a well. Tenv significance labels describe a drug-parent–dose group rather than an individual sample.

**Transcriptome autoencoder**:
A Step 07 model that encodes a measured expression profile and reconstructs that same profile on the frozen CM gene axis.

**Transcriptome representation matrix (转录组表示矩阵)**:
A representation of one measured expression profile with one token per gene in frozen CM gene order. Each token combines the gene's position and its measured expression.

**Shared expression representation interface (共享表达表示接口)**:
The representation generator and restorer used by the Step 07 matrix reconstruction models. The generator maps a profile to its transcriptome representation matrix, and the restorer maps a reconstructed matrix back to expression on the same gene axis.

**Matrix latent representation**:
A fixed number of latent tokens, each with the same number of channels, for one measured expression profile. Its scalar budget is the token count multiplied by the channel count.

**MLPv**:
The Step 07 multilayer perceptron autoencoder with a vector latent representation. In `MLPv_512`, the vector has 512 scalars.

**MLPm**:
The Step 07 multilayer perceptron autoencoder with a matrix latent representation. In `MLPm_512`, the matrix has 512 tokens and 128 channels per token, matching the standard AAE reconstruction latent shape.

**Expression reconstruction tokens**:
The part of the Step 07 matrix latent representation available to the expression decoder.

**Nuisance tokens**:
Extra encoder tokens that can interact with the expression reconstruction tokens but are withheld from the expression decoder. Their role in capturing biological or technical variation remains an experimental hypothesis.

**Training-gene mean** (`mean`):
The Step 07 reconstruction reference that returns one mean expression value per gene estimated from training profiles alone.
_Avoid_: mean Δ, DMSO mean

**Latent XPert** (`latent_XPert`):
The XPert-derived predictor considered for Step 08. It transforms a baseline transcriptome in the latent space of a separately trained transcriptome autoencoder and predicts the drug response on the full CM gene axis.

**XPert landmark reference** (`XPert_base`):
The Step 09 perturbation predictor on the 978 L1000 landmark genes that anchors the gene-space and compression comparisons.

**XPert shared scoring panel**:
The intersection of the genes predicted by all planned Step 09 arms. It
supplies the same tenv drug–dose responses for paired continuous metrics.

**XPert AUPRC@full**:
Direction-agnostic significant-gene retrieval against all significant genes
on the frozen CM axis. Recall includes significant genes outside the
model's coverage in its denominator.

**XPert AUPRC@G**:
Direction-agnostic significant-gene retrieval within a model's prediction
axis. Recall is relative to the significant genes covered by that model.

**Tenv model boundary**:
The parameters updateable together in one backward pass for the tenv prediction
task. A drug encoder trained independently on non-perturbation data is a separate
model and is outside this comparison.

**Shared chemical representation**:
The same structure-derived feature vector supplied for a compound to every
tenv entrant before any model-specific dimension adaptation. Its common origin
does not imply that separately trained adapters produce identical vectors.

**Tenv adapted predictor**:
A perturbation predictor together with its jointly trained linear adapter from
the shared Morgan fingerprint to the predictor's drug-representation interface.
The adapter's parameters belong to the compared model.

**Tenv benchmark sample**:
One model's training-and-test result under one fixed fold and one fixed training
seed in the tenv comparison. The rebuilt comparison takes two folds and two seeds
from the accepted index tables, yielding four run samples per model; the retired
design crossed three folds with three seeds. Neither set is a set of independent
biological samples.
_Avoid_: well, drug–dose group, independent replicate

**Standalone training study**:
A study that trains a model from random initialization and uses a single
source dataset for training, validation and held-out testing. Its validation
and test results describe performance within that source dataset.

**Mix**:
The Step 05 PerturbDiff comparison arm that pools eligible Tahoe and tenv
training data in one fit without balancing the sources.
_Avoid_: Tahoe+CM

**PerturbDiff v2 protocol**:
The uniform nine-slot tenv recipe `perturbdiff-v2-all-validation-ema` defined in
ADR 0042: all-138-group deterministic EMA validation and checkpoint eligibility
every 1,000 steps plus maximum step, with predeclared validation early stopping.
It is distinct from legacy v1's fixed first-40-group validation population.
_Avoid_: exact native recipe, resumed v1

**CMonge v2 protocol**:
The uniform nine-slot tenv recipe `cmonge-v2-full-validation-fixed-draws`
defined in ADR 0042: fixed-draw validation over all 138 conditions with equal
condition weighting and predeclared validation early stopping. Its validation
population differs from legacy v1's single sampled condition.
_Avoid_: exact native recipe, resumed v1

**Joint tenv v2 migration**:
Replacement of all nine CMonge and all nine PerturbDiff formal identities as one
reviewed transition, preserving the other 45 identities and separate legacy
evidence. The formal matrix still expects 63 outcomes. The protected transition
was applied and independently verified: 18 changed identities, 45 unchanged,
22 verified archive operations and no reset of the cumulative allowance.
Application alone does not establish runtime activation or v2 completion.

**Cumulative fold/seed training allowance**:
The original 14,400-second cap shared by all protocol versions and retries of
one tenv model/fold/seed slot. Conservative legacy full-wall debits survive
migration and seed restart; exhaustion yields terminal timeout.
_Avoid_: fresh v2 budget, per-attempt allowance

**Legacy tenv v1 outcome**:
A preserved pre-v2 binary/history/status with its original identity and hashes.
It supplies audit and budget evidence, and cannot fill a formal v2 outcome slot
or establish exact training continuation.

**Tenv significance truth set**:
Genes significant in the measured treatment-versus-matched-control drug–dose
contrast on the full CM axis, retaining their measured up- or downregulation.
This is a statistical significance set rather than a fixed-size response rank.
The rebuilt comparison applies a further declared fold-change threshold; see
*Gated truth gene*.

**Cold-drug fold**:
One Stage 2 benchmark split whose test drugs are absent from its training and
validation sets at every dose and repeat; each drug is tested exactly once
across the complete outer partition. A complete cold-drug comparison covers every
eligible drug exactly once in its test folds; selecting a subset of another
partition does not constitute a complete comparison.
_Avoid_: cold-start split, holdout (unqualified)

**CM-unseen drug (CM 未见药物)**:
A drug whose cardiomyocyte perturbation observations are withheld from model fitting across all CM sources used for adaptation. The drug may have been observed in non-CM pretraining, including Tahoe; this is the prediction setting for drugs without available CM wet-lab responses.
_Avoid_: globally unseen drug, unseen drug (without specifying the CM scope)

**Step 05 cold-cm-drug split**:
A tenv cold-drug split that allows the held-out drug blocks to appear in
Tahoe.

**Cold-cm-fold**:
The Step 09 shorthand for a cold-cm-drug fold: CM test drugs are excluded
from tenv training and validation, while Tahoe may contain those drugs.

**Step 09 cold-cm-drug fold**:
One of the frozen tenv cold-drug folds paired with the shared Tahoe training
and validation partitions. Tenv test drug blocks can appear in Tahoe.

**Stage 2 direction-matched hit**:
A predicted up- or downregulated CM gene that is significant in the measured
drug–dose contrast and has the same direction as its measured log2 fold change.

**Stage 2 AUPRC**:
The macro mean of exact group average precision for combined direction-matched
retrieval, with every measured significant gene on the 18,164-gene CM axis in
the recall denominator, including genes outside prediction coverage.

**Stage 2 directional AUPRC**:
`AUPRC_positive` retrieves significant upregulated CM genes using positive
predictions; `AUPRC_negative` retrieves significant downregulated CM genes using
negative predictions. Each uses its full-CM directional truth denominator.

**Stage 2 pooled PR curve**:
Precision and recall of pooled gene–group predictions across a fixed eligible
group population at common prediction-magnitude thresholds. Its endpoint is
the attainable pooled recall; its area differs from macro AUPRC.

**Step 05 cold-drug split**:
A tenv cold-drug split that also withholds its tenv validation and test drug
blocks from Tahoe training and validation.

**Step 05 warm-start split**:
A tenv split that holds out complete parent–dose combinations and retains
other doses of each held-out parent in training.

**Δ space**:
The CM benchmark target: `log1p_cp10k(treated) − mean(log1p_cp10k(same-plate DMSO))`,
on the ordered 18,164-gene CM axis. It expresses a treated well's change from its
matched control reference; group targets average these well-level changes.
_Avoid_: fold change, expression delta (unqualified)

**Group**:
One `(drug_key, dose)` pair in a CM benchmark, represented by the mean response
of its measured wells. Groups receive equal weight in benchmark metric macro
averages.
_Avoid_: condition (unqualified), sample

**Paired Transfer comparison**:
A Step 06 comparison between original Tahoe-to-tenv Transfer and one
intervention on the same tenv drug–dose groups. Each arm is evaluated at
its own author-loss-selected EMA checkpoint.
_Avoid_: combined optimization arm

**Step 07 transfer pair**:
Two results for the same model and seed on the shared tenv fold, one from
the Tahoe-pretrained predictor and one from its tenv-fine-tuned descendant.

**Cold-drug-capable model**:
A model whose drug input is derived from chemical structure (fingerprint, RDKit descriptor, or a structure encoder), so it can score a drug it never saw in training. Stage 2's admission gate.
_Avoid_: warm-drug model, structure model

**Entrant set**:
The models admitted under a particular benchmark's frozen task, code and
resource gates. The historical seven-model set does not establish admission to
the redesigned tenv comparison.
_Avoid_: candidate list, all models

**Ideal-fit label**:
The `ideal_fit` report column: `full` (author architecture is used as intended under the contract), `partial` (a declared omission or transductive control, e.g. chemCPA's skipped 978 pretrain, biolord's transductive prediction, PerturbNet's missing control transcriptome), or `none` (a CPU baseline). One leaderboard; the label is not a second table.
_Avoid_: ideal model, tier

**Benchmark run**:
One `(model, fold)` unit in the historical ADR 0027 benchmark: fit on the training drugs, select the checkpoint on the cold validation drugs, predict the test wells, and write the per-unit artifacts and `status.json`.
_Avoid_: experiment, trial


### Stage 2 tenv model comparison

**Tenv dose-input comparison**:
The Stage 2 comparison of ways to supply continuous dose to TranSiGen across
the measured tenv doses, with unseen-drug prediction as its primary objective
and held-dose prediction as an auxiliary objective.

**Dose-input round1**:
The Stage 2 TranSiGen dose-input optimization exploration. Its ADR 0078
successor compares seven fixed configurations independently on Tahoe and tenv
using axis-local G_5000, with cold-drug and cold-drug-dose evaluation. The
historical ADR 0076 generation used tenv with cold-drug and held-dose evaluation.
Both generations provide comparison evidence without a downstream weight release.
_Avoid_: final deployment selection, universally optimal dose architecture

**Dose-input comparison run**:
One source, fixed method, scenario, fold and seed in the dose-input comparison.
A historical held-dose fold denotes a distinct omitted dose. Runs do not
represent independent biological replicates.
_Avoid_: independent replicate, Benchmark run (the historical ADR 0027 unit)

**Cold-drug-dose**:
A dose-input round1 test at the intersection of unseen drug parents and one
unseen measured dose. The development partition contains other parents at
the other doses. The remaining two parent-by-dose strips are excluded.

**Sample-equal Delta evaluation**:
Calculate each treated sample's gene-wise Delta metrics with its own matched
control, then give each sample equal weight within a reported population.
This is the ADR 0078 evaluation unit for both sources. The historical ADR 0076
comparison retains parent-equal evaluation.

**Parent–dose response**:
The mean transcriptional response of repeated treated wells for one drug
parent at one measured dose, retaining that dose as a distinct evaluation unit.

**Tenv comparison step** (`02_models_tenv-5k`):
The active Stage 2 comparison of perturbation predictors on single-dose `tenv`, with parent-level held-out Δ performance as the primary result. It owns the shared comparison splits, entrants, predictions and figures; the released inputs and palette stay in `01_bases`.
_Avoid_: benchmark step (when meaning this comparison)

**Direct metric**:
`PCC`, `Spearman`, `RMSE` or `R²` comparing predicted and observed Δ across the gene axis. The current primary unit is a drug parent after out-of-fold repeated-well aggregation; per-well metrics are auxiliary.

### Historical Stage 2 reliability vocabulary

The following terms describe the retired ADR 0064 comparison and its sealed records.

**@d reliability family** (reliability family):
`Precision@d`, `Recall@d`, `F1@d`, `Spearman@d` and `Reliability@d`, computed per drug–dose group from that group's mean predicted Δ against its measured DESeq2 contrast. It replaces the retired `k`-parameterized indirect family. The `@d` suffix is the family marker: every member is parameterized by the predicted-effect threshold `d`, not by a declared table size.
_Avoid_: indirect metric (unqualified, when the @d family is meant); batched metric (withdrawn, never built)

**Measured effect**:
A gene of one drug–dose group whose measured contrast has `padj < 0.05`; its direction is `sign(log2FC)`. It is the label side of the @d family.
_Avoid_: gated truth gene (the retired `|log2FC| >= 2` gate), DEG (unqualified)

**Predicted effect threshold** (`d`, `d'`):
The threshold applied to a group's `|Δ̂_mean|`: a gene is *predicted changed* when `|Δ̂_mean| >= d`. On a run it is declared as `d = k · 0.01 · M` for `k = 1..100`, where `M` is that run's largest `|Δ̂_mean|` over its test groups and the whole gene axis; `d' = d/M = k/100` is the plotted, scale-free axis. `d` is a predicted effect magnitude, not a fold change and not a p-value.
_Avoid_: cutoff, significance threshold (a predicted magnitude carries no statistical test)

**Threshold window** (`W(d)`):
The genes of one group whose `|Δ̂_mean| >= d`. It is the model's own "changed gene" table at that threshold, so its size follows the group's predicted effect strength instead of a declared table size.
_Avoid_: predicted set, top-k

**Direction-correct hit** (`H(d)`):
A gene of `W(d)` that is a measured effect and whose measured direction equals the predicted direction `sign(Δ̂_mean)`.

**Precision@d / Recall@d / F1@d**:
`Precision@d = |H|/|W|`, `Recall@d = |H|/|measured effects|`, and `F1@d` their harmonic mean. All three are total by convention: an empty window on a group with no measured effect scores `F1 = 1`, an empty window on a group with measured effects scores `0`, and a non-empty window on a group with no measured effect scores `0`. Harmonic rather than arithmetic averaging is load-bearing: an arithmetic mean peaks at the smallest table and would make `Reliability@d` select a one-gene table.

**Spearman@d**:
The Spearman correlation inside `H(d)` between predicted magnitude (`|Δ̂_mean|` rank) and measured magnitude (`|log2FC|` rank): given that those genes are right and correctly directed, did the model order them by effect size correctly? Neutral (`mult = 1`) when `|H| < 3`, and also when the predicted or measured magnitudes of `H` have no rank variation; both rules can only inflate a score, so they are declared here and in ADR 0064. It is not the retired `Spearman@k`, which correlated signed predictions against measured Δ over all of a top-k table.
_Avoid_: Spearman@k (the retired top-k statistic)

**Reliability@d**:
`F1@d · mult`, with `mult = (1 + Spearman@d)/2` for `|H| >= 3` and `1` otherwise. The macro `Reliability@d` of a run is the equal-weight mean of the per-group values, **not** `F1@d · (1 + Spearman@d)/2` applied to the macro `F1@d` and macro `Spearman@d` columns; the two agree only when Spearman is constant across groups. Range `[0, 1]`: `0` means the emitted table is wrong or direction-inverted, `1` means it is exactly the measured effect set in the correct order. It is a model-selection proxy for the downstream phenotype task, not a biological claim.

**`d*`**:
The `d` that maximizes a run's macro `Reliability@d` over the 100-point scan. There is one `d*` per run, so a consumer of a given run's checkpoint uses that run's `d*`.

**Predicted absolute rank table (预测绝对秩表)**:
One group's genes listed as `gene_index`, `rank`, `direction` and `|Δ̂_mean|`, ordered by predicted absolute effect. A light rank table orders the same genes by absolute *observed* Δ.
_Avoid_: prediction table, prediction matrix

**Mean baseline** (`Mean`):
The reference entry that predicts, for every gene of a group, the mean Δ of that same group in the training partition. It is deterministic, not a trained fit, and its fold copies are not independent training results.
_Avoid_: delta_mean (retired name), no_change (withdrawn)

**Fold/seed index table**:
`01_bases/output/fold_n_seed/folds.csv` or `seeds.csv`. Each table orders frozen identities for studies that explicitly consume that table. A separately frozen comparison partition has its own fold identities; taking a prefix of a different partition does not recreate it.

**Stage palette**:
The two-dimensional `01_bases/output/palette/` table of plotting colors: six hue families (gray, blue, green, orange, purple, magenta) crossed with nine lightness-and-saturation shades. A color means the same thing in every Stage 2 figure.

**Learning-curve figure**:
`figure/learning_curve.png`: one subplot per learnable model, epoch on the horizontal axis and loss on the vertical axis. The blue family is training and the green family is validation; one shade identifies one `(fold, seed)` run in both families.

**Reliability figure** (`figure/metrics_indirect.png`):
Two rows by four columns: `Precision@d`, `Recall@d`, `Spearman@d` and `Reliability@d` as curves over `d'` (with `F1@d` as a thin reference line on the Reliability panel), then `Reliability@d*` and `d*` bars, a macro three-class confusion heat map at `d*`, and the null-versus-changed Reliability split.


### Stage 2 gene axis

Vocabulary for the current ADR 0073 full-source tier-first DEG1000/f1000 release. Earlier ADR 0071/0072 measured contracts are historical.

**Gene axis** (`G_N`):
The ordered gene set a Stage 2 perturbational model reads and predicts, as an ascending
subset of the frozen 18,164 gene indices without renumbering. ADR 0073 specifies sizes
2,000, 3,000, 4,000 and 5,000, all present in the measured release.
_Avoid_: narrowed corpus axis, new identity axis

**Gene-axis candidate**:
A frozen-axis gene under consideration for the axis — the rows of
`output/gene-axis/candidate_rank_sum.csv`, one per frozen gene. tenv-HVG ranks all 18,164
genes, so every frozen gene is a candidate; a candidate is not thereby a DEG. Under ADR 0073,
only membership in an eligible disease's DEG1000 gives disease-side support; HVG allocation
and fallback can still admit other candidates.
_Avoid_: gene universe, axis pool

**{disease}-DEG**:
A Control-vs-Disease differential-expression result from a source admitted under ADR 0071's
tiers — CM pseudobulk from single-cell data, CM-nuclear bulk, in-vitro CM/organoid/engineered
tissue, or tier-3 human atria/ventricle bulk
tissue that is not cardiomyocyte-specific; ADR 0073 removes the former cumulative-count
condition on collecting lower tiers. A whole-tissue bulk series is a tier-3 DEG source,
never a tier-1 or tier-2 one; ADR 0071 supersedes ADR 0066's blanket exclusion of it.
_Avoid_: disease source (unqualified, when a particular tier is meant)

**TOX** (ADR 0072, retained by ADR 0073):
The single small-molecule drug-cardiotoxicity disease category: transcriptional responses to drugs with independently supported cardiotoxicity, across qualifying exposure conditions, regardless of mechanism class. Positivity belongs to the drug and does not assert toxicity at every dose or time.
_Avoid_: cardiotox_TOP2 / cardiotox_other (as successor disease categories), all drug exposures, toxicity-positive exposure (when drug-level admission is meant)

**Disease eligibility** (ADR 0072, retained by ADR 0073):
A declared disease's qualification to contribute to cross-disease fusion: its distinct DEG union on the frozen CM gene axis contains at least 100 genes. An excluded disease still has its source evidence, DEG analysis, complete union, fusion and DEG1000 record.
_Avoid_: source eligibility (a different decision), disease deletion, AF-specific exclusion

**Disease category**:
A grouping whose datasets are analyzed and fused together. ADR 0072 declares
`nHF`, `ICM`, `AF`, `genetic_cardiomyopathy` and `TOX`; eligibility is a separate
decision. The measured ADR 0073 release has four eligible terms. AF retains its 25-gene
DEG1000 and upstream evidence but contributes no cross-disease term.
The pre-successor ADR 0071 release had six terms including two mechanism classes.
_Avoid_: disease label, cohort, source (a category is not a source)

**nHF**:
The non-specific human heart failure category: qualifying clinical non-ischemic CM evidence, compatible non-genetic CM models and human ventricular bulk under the source-tier rules. Explicit genotype-split cohorts include the Reichart donors without a declared pathogenic genotype; other source-specific genotype uncertainty remains declared rather than inferred.
_Avoid_: NICM (when the category is meant)

**ICM**:
The ischemic cardiomyopathy category, merging chronic ischemic heart failure and acute MI
samples into one ischemic continuum because no human deposit verifies a per-sample acute
window.
_Avoid_: IHF, MI/MIRI (the rejected chronic/acute split)

**genetic_cardiomyopathy**:
The hereditary myocardial disease category: genotype-positive human donors plus in-vitro
isogenic hiPSC-CM models of HCM, DCM (TTN/LMNA/RBM20), ARVC and DMD.
_Avoid_: hereditary CM (unqualified)

**cardiotox_TOP2** / **cardiotox_other**:
Historical pre-successor drug-cardiotoxicity categories (ADR 0066 and ADR 0071): topoisomerase-II poisons, and every other mechanism
class holding at least three small-molecule drugs. Antibody and peptide drugs are excluded.
_Avoid_: chemotherapy toxicity (unqualified)

**Source rank**:
A dataset's within-dataset ordering of its **DEGs** by the effect score
`|log2FC| × −log10(padj)` — DEGs rank `1..n_DEG`, every other frozen-axis gene gets
`n_DEG + 1`, and `rank_norm = rank / n_DEG`. The denominator is the DEG count, never the
measured-gene count: a dataset with three DEGs supplies three ordered genes and a wall of
`(n_DEG + 1) / n_DEG`, not a complete 18,164-gene ordering (ADR 0071 decision 3, superseding
ADR 0066/0068's measured-gene rank). A source with `n_DEG = 0` has no `rank_norm` at all and
is not a fusion source. tenv-HVG is not a DEG source and is never ranked this way — it
contributes its dispersion rank to the axis's first half only. Ranks are never compared
across datasets as if they were effect sizes.
_Avoid_: global rank, meta-rank, measured-gene rank (retired)

**Class source**:
The single disease term derived from its informative source ranks. ADR 0073 orders the complete DEG union by highest DEG source tier, within-disease rank sum, then gene index, before retaining DEG1000.
_Avoid_: merged dataset, meta-analysis

**Rank-sum fusion**:
The cross-disease score `rank_sum = Σ_d rank_d / n_d` over eligible diseases, where `n_d` is the complete axis-resolved disease DEG union size. Under ADR 0073, a gene inside DEG1000 carries its tier-first disease rank and any outside gene carries `n_d + 1`; tenv-HVG is not a term.
_Avoid_: weighted fusion, truncated-set denominator, seven-source sum

**pN** (retired 2026-11-15):
The v1 cross-source replication count: how many of the fused sources (five categories plus
tenv-HVG) had a gene in their top-`N` set, for `N ∈ {2000, 5000}`. It is **gone** — the
`p2000`/`p5000` columns, their readers, their tests and the measured-gene rank they were
computed from are removed (ADR 0071 decision 9). Its replacement was `f2000` in ADR 0071/0072; ADR 0073 replaces that with `f1000`. Both count only
eligible diseases and never count tenv-HVG. Any v1 text,
figure or log that still shows `pN` is a sealed historical record; do not reintroduce the
term into a live product.
_Avoid_: pN (for the live signal), support count, vote

**Axis fill rule** (ADR 0073):
For N in {2000, 3000, 4000, 5000}, select the top N/2 tenv-HVG genes, then N/2 unused candidates with f1000 at least 1 ordered by f1000 descending, cross-disease rank sum ascending and gene index ascending; fill any shortfall with unused HVG genes in HVG order. Publish each exact-size axis in ascending frozen gene index; the four axes are nested.
_Avoid_: DEG-only whole axis, non-nested axes, top-N selection (unqualified)

**Source top-N coverage**:
For tenv-HVG, the fraction of its top-N genes present in the N-gene axis; for an eligible disease under ADR 0073, the fraction of its own DEG1000 present in that axis. Excluded diseases have unavailable coverage with an explicit exclusion reason, rather than zero coverage; historical releases retain their DEG2000 denominators.
_Avoid_: source recall, contribution quota, shared disease denominator

**Label-permutation null**:
The null distribution of a dataset's axis-ordering statistic obtained by permuting its sample
group labels (100 permutations, same DE model refit, top-`N` overlap with the released axes); a
dataset newly admitted to the fusion must beat it. A uniform-random gene baseline is not a valid
substitute: it flatters any dataset whose measurable genes happen to be large-effect genes (ADR
0070).
_Avoid_: random-gene null, permutation p (unqualified)

**Source dataset table**:
The category subtables at the start of `stage2/01_bases/gene_axis_fill_rule.docx`: one row per released dataset record with its
category, accession, sample type, groups, sample sizes, sequencing technology and the published
v2 DEG count (`catalog.json` `n_deg`, never the padj-only `n_significant`), plus a separately labelled tenv-HVG non-disease row. Shared studies/controls and sample units remain explicit.
_Avoid_: dataset summary, metadata table

**Inflection point**:
The candidate `N` marked by the rank-sum curve's Kneedle plus second-difference criteria,
reported with its parameters and stability and allowed to be `"inflection": false`. It is
diagnostic only: the axis sizes are declared, not read off the curve.
_Avoid_: optimal N, elbow (unqualified), cutoff

**Gene-axis view**:
A View whose expression columns follow an ordered subset of original frozen gene indices; retaining a gene-axis list does not imply publishing its expression view. ADR 0077 retains four named axes but publishes expression only on G_5000.
_Avoid_: corpus (for a view), subset matrix

**Axis-local normalization**:
Per-profile CP10K followed by natural log1p, with the library total calculated only from raw counts on that profile's declared gene axis. Matched-control expression is the mean of individually normalized vehicle profiles; full-axis normalization instead uses the complete frozen gene axis as its denominator.
_Avoid_: re-normalizing log expression, full-axis normalization (for this denominator)

**Raw-count release** (`raw_counts/`):
tenv's native-gene raw counts with probe evidence and Tahoe's token evidence — the replayable
bridge from immutable `source/` bytes to the frozen axis.
_Avoid_: raw corpus

**{N}k eval truth**:
The evaluator-only ground truth a view will own (`2k_eval_truth/`, `5k_eval_truth/`) once the
model-evaluation design is reworked; no truth is released until then.
_Avoid_: tenv-truth (the retired full-axis release)


### Stage 1 evaluation

**Top-k true set**:
The ground-truth label convention of the Stage 1 XPert benchmark — the step-05 official folds and the `warmstart_pr/` PR/AUPRC analysis: per sample, the true set is the top-k up genes and the bottom-k down genes of that sample's true Δ (k = 100 in the 978-landmark space), every other gene being a negative in both directions. A fixed-truncation (ranking) definition with no significance test, so each direction's prevalence is exactly k/978 and that prevalence is the random AUPRC baseline.
_Avoid_: DEG set (unqualified), differential set

Historical Stage 2 v1 benchmarking terms (direction-aware AUPRC, unified test set, data-ablation benchmark, dtoxs-primary benchmark) are retained in Git history with the former `stage2/_archive/cm_bench_v1/` and are not active vocabulary.

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
