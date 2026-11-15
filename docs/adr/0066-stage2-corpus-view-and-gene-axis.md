# Stage 2 01_bases: raw counts, corpus, gene-axis views and the two-version gene axis

**Status:** accepted (2026-11-13). The user re-scoped the Stage 2 foundation and fixed the
gene-axis design in the 2026-11-13 design session. This record supersedes ADR 0065 in full,
supersedes the single-corpus `X`/`control`/`delta` semantics of ADR 0061 and the
`01_bases/output/` layout contract of ADR 0062 for this step, and retires the measured-truth
and successor-envelope products together with this step's committed `log/` convention. ADR
0030's external-item provenance contract and ADR 0062's candidate-then-replace acceptance
rule still apply.

*Amended later, text unchanged (accepted records are historical).* ADR 0068 changed the
fusion's non-informative-source rule and rank normalization; ADR 0069/0070 removed the `AF`
category; **ADR 0071 supersedes this record's *Source eligibility*, *Priority* and
*Candidate table and the two-level rank-sum fusion*** — the measured-gene rank, `pN`, the
seven-fusion-source framing (six categories plus tenv-HVG) and the exclusion of whole-tissue
bulk — and the corpus/view semantics above it remain in force.

## Why this record exists

Two things changed at once.

1. **The "HVG" concept in this step was ambiguous and partly wrong.** One term covered two
   different products computed on two different populations, and the disease side of the
   gene-axis evidence admitted whole-heart bulk tissue that is not cardiomyocyte-derived.
   The vocabulary is corrected below, and the correct vocabulary is now the contract.
2. **The axis is no longer undecided.** Stage 2's exit is the predicted Δ on a *selected*
   gene axis (`docs/roadmap.md`); ADR 0065 deliberately prepared evidence without choosing
   `N`. The user has now fixed `N ∈ {2000, 5000}` and a 50/50 selection quota, so this step
   releases **two axis versions**, and it rebuilds the `output/` tree so that raw counts,
   the model-facing matrices and the axis evidence stop sharing one directory.

ADR 0065's "no axis release" boundary, its unweighted quota-free fusion, its direction
composition, its tenv non-HVG curve and its animal half-weight column all disappear with it.
The measured truth and the successor envelope disappear too: the model-evaluation design is
being reworked, and each view will build its own evaluator truth once that design exists.

## Vocabulary: what "HVG" means here

The two sources of the gene axis are now named for what they are, and neither name is reused
for anything else. Both are resolved to the frozen 18,164 latest-HGNC symbols, so the two
sources and the released axes share one symbol space.

**tenv-HVG** — the **per-well expression dispersion of the tenv corpus**. It is computed over
the **4,206 treated wells** of `tenv` (every chemical parent × dose × replicate; vehicle
controls are not rows), from the frozen CP10K/log1p expression, as the published normalized
dispersion: sample variance/mean of linear CP10K recovered by `expm1`, log dispersion
z-scored in 20 equal-width `log1p(mean)` bins, singleton bins scored 1, undefined sorts last,
`gene_index` tie-break. The published product is `output/tenv-hvg/descriptive.csv`. It
measures how a cardiomyocyte transcriptome moves **in general** under chemical perturbation,
which is why it is useful for suppressing training-data noise; it carries **no direction**.

An intermediate proposal to extend this population to all 4,478 tenv wells (adding the 272
vehicle wells) was **rejected**: the axis must not create a third tenv population, and the
published ranking is already accepted and replayable.

**tenv-HVG (model selection)** — the train-only fold rankings in `output/tenv-hvg/`
(`fold_0/1/2.csv`). They are a *different population* used by model gene selection under ADR
0061 and are **not** an axis input. Anything that says "tenv-HVG" without a qualifier means
the axis population above; model selection must say "tenv-HVG (model selection)".

**{disease}-DEG** — a Control-vs-Disease differential-expression result from a
**cardiomyocyte-derived** source. It is the disease input of the axis and is produced per
dataset under the eligibility and class rules below.

## Source eligibility (this is the correction)

A {disease}-DEG source must be **cardiomyocyte-derived**, in one of exactly two forms:

- a **bulk** DEG from in-vitro cardiomyocytes, cardiac organoids or engineered cardiac
  tissue; or
- a **cardiomyocyte pseudobulk** DEG from single-cell data, i.e. the CM cells of a donor are
  summed into one profile before differential testing.

Whole-heart or whole-left-ventricle **tissue bulk** is **not** eligible, even though it
contains cardiomyocytes: its signal is dominated by cell-composition and fibrosis shifts that
do not describe the cardiomyocyte. This disqualifies the whole-tissue sources the previous
evidence round used — `GSE141910`, `GSE116250`, `GSE128188`, `GSE282618` (human LV/atrial bulk)
and `GSE160516`, `GSE108940` (mouse whole-LV bulk MIRI).

**Priority.** For each declared disease category, the highest available tier supplies the
source:

1. human-sample single-cell (CM pseudobulk);
2. hiPSC/hESC-CM, and organoids or engineered tissue built from them (bulk);
3. animal-sample single-cell.

All independent datasets **at the winning tier** are kept; they are independent cohorts, not
duplicates, and their agreement is what the fusion is for. Shares controls are declared (see
below), never pooled as extra replicates.

**Isogenic exception.** In-vitro models carrying a patient genotype with an isogenic control
are kept as their own evidence type even when a human single-cell cohort exists for the same
category: a mixed-etiology clinical cohort and a single-genotype isogenic model answer
different questions, and collapsing them loses the causal signature.

**Drug cardiotoxicity.** The unit is a **mechanism class**, not a drug. A class enters only
when the source contains **at least three small-molecule drugs** of that mechanism; antibodies
and peptides are excluded (`trastuzumab`, `bevacizumab`, `rituximab` drop out). Class members
are pooled into one Disease group and compared with the same-plate vehicle.

**Genotype.** Where a source declares per-donor genotype, the donors are split by it rather
than by the study's disease label.

## The category roster

Six disease categories plus the tenv-HVG source give **seven fusion sources**. A category is
one source; its datasets are fused inside it (see next section).

| category | source datasets (all at the winning tier) |
|---|---|
| `nHF` — non-specific human heart failure: clinical, non-ischemic, non-hereditary | `GSE270788` NICM arm (7 donors vs 6 Donor); the Reichart cellxgene CM object's donors without a declared pathogenic genotype |
| `ICM` — ischemic cardiomyopathy, chronic and acute merged | `GSE270788` ICM arm (5) and `GSE270788` AMI arm (5), each against the same 6 Donor controls |
| `AF` — atrial fibrillation | `GSE255612` (18 AF / 16 control donors, left-atrial CM) |
| `genetic_cardiomyopathy` — hereditary myocardial disease | Reichart genotype-positive DCM donors and its 8 ARVC donors (human CM pseudobulk; its single NCCM donor is not used); isogenic hiPSC-CM bulk models `GSE224129` (HCM, MYBPC3, EHT), `GSE234141` (DCM, TTN truncating), `GSE304575` (DCM, LMNA H222P), `GSE305648` (DCM, RBM20, cardioid), `GSE253223` (ARVC, PKP2 knockdown), `GSE194297` (DMD, cardiac organoid) |
| `cardiotox_TOP2` — topoisomerase-II poison cardiotoxicity | `GSE243674` (DNR/DOX/EPI/MTX), `DToxS-iPSC`/GSE217421 (DAU/DOX/EPI/IDA), `GSE244740` (DOX/IDA/Mitox) |
| `cardiotox_other` — the remaining cardiotoxic mechanisms | `DToxS-iPSC` non-TOP2 mechanism classes that hold ≥3 small-molecule drugs, each class its own DEG inside the category |

**Why ICM merges chronic and acute.** An earlier split into `IHF` (ischemic HF, non-acute,
> 2 weeks after the insult) and `MI/MIRI` (acute, ≤ 2 weeks) was **rejected**: no human
single-cell deposit carries a per-sample time-since-infarct field that can be verified
against a 14-day boundary, and the animal single-cell survey found no chronic ischemic series
with CM labels, a sham arm and ≥ 3 replicates. A category boundary that cannot be measured is
not a category boundary. `GSE270788`'s AMI arm enters `ICM` as the acute end of one ischemic
continuum, and its author-declared "acutely infarcted" status is recorded as
`time_window_basis: author_declared`, not as a verified window.

**Why nHF is "non-specific".** `GSE270788`'s NICM arm carries no genotype column, so
"non-hereditary" cannot be asserted donor by donor; it is a clinical, non-ischemic cohort and
is recorded with `genotype_status: not_declared`. The Reichart object does carry
`Primary.Genetic.Diagnosis`, so there the genotype split is applied and the positive donors
move to `genetic_cardiomyopathy`.

**Declared gaps.** `HHD` (no cardiomyocyte-derived human single-cell candidate), `PPCM` (no
human single-cell hits), ARVC patient missense mutations, HCM `MYH7`, LQT/CPVT, and any
animal single-cell replication. Gaps are reported, not filled with ineligible tissue data.

**Sources deliberately not used.** `disease_data/GSE137293` (pig scar-vs-border, a severity
gradient), `pert_data/LINCS-LDS-1593` (L1000-978 coverage note), the myocarditis cellxgene
object (it shares its control donors with the Reichart object, so the two cannot both count as
independent), `GSE183852` and `GSE298023`/`GSE145154` (available replication, but cell labels
live only inside multi-GB `Robj`/`rds` files), and `pert_data/GSE247541` (recorded outside the
Stage 2 union scope by its item policy; admitting it needs its own decision).

**Download budget.** Nine new external items (the three human single-cell deposits, six
isogenic models), about 5.1 GB; the three cardiotoxicity sources are already owned read-only (`DToxS-iPSC` and
`GSE244740` under `pert_data/`, `GSE243674` under `disease_data/`). This is inside the declared ceiling of 12 items and 100 GB, all open access.

**Shared-control declaration.** All three `GSE270788` arms use the same six Donor controls;
the `DToxS-iPSC` classes share the same vehicle wells; Reichart's controls are shared with
the dropped myocarditis object. Each source records a `shared_control_group` identifier. The
fusion does not correct for this and claims no independence it does not have.

## Candidate table and the two-level rank-sum fusion

The candidate pool is the frozen 18,164-gene axis (ADR 0061); no gene outside it can enter,
because no model reads outside it. tenv-HVG ranks all 18,164 genes, so the pool is complete.

**Per-dataset DEG rank.** Within one dataset, and within that dataset only, a gene's rank is
its position in the descending continuous effect score `|log2FC| × −log10(padj)` over the
genes the table measures; a multi-contrast source collapses per gene by the maximum effect;
`NA` `padj` counts as non-significant; a gene the dataset does not measure gets
`n_measured + 1`.

**Rank normalization.** `rank_norm = rank / n_measured`, where `n_measured` is the number of
genes the table measures — **not** the number of significant DEGs. Normalizing by the
significant count would inflate a dataset that found almost nothing to the same weight as a
dataset that moved half its transcriptome.

**Level 1 — inside a category.** Every dataset at the winning tier produces its own DEG; the
category's rank is the ascending order of `Σ rank_norm` over its datasets, ties broken by the
frozen `gene_index`. A single-dataset category is simply its own order.

**Level 2 — across categories.** `rank_sum = Σ rank_norm` over the seven sources (six
categories + tenv-HVG), ascending, ties broken by `gene_index`. The **master table**
`output/gene-axis/candidate_rank_sum.csv` carries `symbol`, `gene_index`, `tenv_hvg_rank`, one
`<category>_rank` per category, `rank_sum`, `n_sources_measured`, `p2000` and `p5000`.

**pN.** `pN` is the number of the seven sources whose top-`N` set contains the gene
(`N ∈ {2000, 5000}`). It measures cross-source replication: how many independent disease or
variability mechanisms place the gene near the top of their own list.

**Inflection is diagnostic only.** The rank-sum curve still gets the declared Kneedle +
second-difference check and may honestly report "no inflection"; it does not choose `N`. The
histogram of `rank_sum` over all candidates, with both axes' intervals marked, is the
evidence panel.

## The two released axes

For each target size `N ∈ {2000, 5000}`, `G_N` is built in two halves:

1. the top **`0.5·N`** genes of tenv-HVG (1,000 for `G_2000`, 2,500 for `G_5000`);
2. then **`0.5·N`** more from the remaining candidates, ordered by **`pN` descending, then
   `rank_sum` ascending, then `gene_index` ascending**. Lower rank sum is better; ordering by
   descending rank sum would select the worst genes, which was an explicit ambiguity resolved
   in the design session.

The two axes are **not nested** — different `pN`, different tenv top sets — and the manifest
reports their intersection so no consumer treats one as an increment of the other.

An axis is an **ordered subset of the frozen `gene_index` space**: it keeps the original
integers and is ordered by ascending `gene_index`, never renumbered to `0..N−1`. Integer gene
identity is the repository's error-detection device (ADR 0061), and a renumbered axis would
silently change what an index means.

## The `01_bases/output/` layout

```
output/
  api/            # the single public reader/metric implementation (ADR 0062)
  palette/        # shared figure palette
  fold_n_seed/    # fold/seed partitions
  gene-identity/  # frozen 18,164 identity tables
  raw_counts/     # native raw counts and mapping evidence (tenv 19,613 native genes;
                  #   Tahoe token evidence) — the `source/` → axis bridge, replayable
  corpus/         # 18,164-axis integer raw counts for tenv and tahoe_bulk, with metadata,
                  #   control correspondence and role flags; no normalization, no Δ
  tenv-hvg/       # published descriptive + train-only rankings
  disease-deg/    # per-dataset DEG tables and their within-dataset ranks (level-1 inputs)
  gene-axis/      # candidate_rank_sum.csv, gene_axis_2000.csv, gene_axis_5000.csv, manifest
  view-2k/        # model-facing 2,000-gene view: X, control, Δ
  view-5k/        # model-facing 5,000-gene view: X, control, Δ
```

`view-2k/`/`view-5k/` are **views owned by this step** and consumed by the model steps; after
the evaluation redesign each view will additionally own its `{N}k_eval_truth/` under `output/`.
That directory does not exist yet and no truth is released by this record.

## Corpus and view semantics

**`corpus/`** holds integer counts on the frozen 18,164 axis. Row sets: `tenv` 4,206 treated +
272 vehicle wells; `tahoe_bulk` 65,318 treated + 1,400 `DMSO_TF` rows. Every row carries
`role ∈ {treated, vehicle}`, the source profile/well identity, plate and cell identity, the
desired parent structure, the micromolar dose and the matched-control group id. Controls are
**rows**, not a separate lookup. The corpus carries no normalization, no Δ, no HVG subsetting
and no filter: those are the view's job, so that one raw truth serves every axis.

**A view** is built from the corpus in three deterministic steps:

1. normalize each profile on the **full 18,164-gene library total** and take `log1p` of CP10K
   (the frozen ADR 0061 normalization, unchanged);
2. subset to the axis's genes;
3. attach `control` = the mean of that row's matched controls normalized the same way, and
   `Δ = X − control`.

**The denominator is the full axis, not the view.** Re-normalizing inside the view is the
obvious alternative and was **rejected on measured evidence**. On tenv, with the 2,000-gene
axis carrying 12.9% and the 5,000-gene axis 36.3% of the full-axis CP10K, view-local
normalization would make the per-gene input baseline differ by 2.5× between the two axes and
`mean|Δ|` differ by 1.87× (0.666 vs 0.355) purely as an artifact of the denominator. Since the
point of releasing two axes is to compare `PCC`, `Spearman` and `RMSE` of the same model across
them, `RMSE` must share one scale; the view-local choice also forces per-axis
hyperparameter retuning and adds a sum-to-constant constraint that biases per-sample
gene-wise correlation. The measured comparison is recorded in the design record; the
decision is: **one normalization, two gene subsets**.

**Tuning protocol for the axis comparison.** Because the two views now share a scale, the
comparison must still declare its protocol: one hyperparameter set held fixed across axes
(the clean axis comparison), or a per-axis tuning budget (each model at its best). Both are
reported, together with the metrics restricted to the shared genes of the two axes, so that a
difference in `RMSE` can be separated from the difference in gene-set composition.

## Removals

Products removed by the rebuild, after the new products are accepted: `shared-expression/`,
`training-corpora/`, `tenv-raw-counts/` and `tahoe-raw-counts/` (their content moves into
`raw_counts/`), `tenv-truth/`, `tenv-truth-successor/`, `historical-training-release/`, and
the previous `gene-axis/` file set (`source_summary.csv`, `fused_rank.csv`, `overlap.csv`,
`inflection.json`, `tenv_non_hvg_curve.csv`, `direction_summary.csv`,
`mouse_human_orthologs.csv`). The mouse ortholog table and its R builder go because no animal
source survives eligibility.

Code removed: the truth and evaluation chain (`build_significance.py`, `significance.py`,
`deseq2_contrasts.R` and its sealed copy, `truth_successor.py`, `bind_truth_successor.py`,
`reissue_truth_successor.py`, `evaluate_predictions.py`, `replay_pr.py`,
`verify_hvg_metrics.py`, `export_diagnostic_predictions.py`, `publish_diagnostic_figure.py`,
`bases/evaluation.py`, `bases/metrics.py`, `bases/plots.py`, the diagnostic PR figure), the
one-off migration and audit scripts (`migrate_products.py`, `audit_metadata_bindings.py`,
`verify_payload_digests.py`), the superseded build chain (`build_shared.py`,
`build_corpora.py`, `verify_training_source.py`) and `build_mouse_orthologs.R`. Git retains
all of it.

**`log/` is deleted, not archived.** In this repository `log/` is execution-run intermediate
output and verification results, not a history of directory changes (git is that history).
It is therefore regenerated at acceptance and never kept across a refactor. The step's
`log/` tree goes with the rebuild.

## Deferred

- `02_models_tenv` reads `training-corpora/tenv.h5ad` layers `control` and `delta` directly;
  it must be re-pointed at `view-2k/`/`view-5k/`, and its accepted 2-fold × 2-seed matrix
  re-run. That redirect and re-run are **not** in this record's scope.
- The evaluation redesign ("模型评价方式") is not in scope; when it lands, each view creates
  its own `{N}k_eval_truth/` under `output/`.
- `GSE183852`, `GSE298023`, `GSE145154`, the myocarditis object and the excluded drug data
  stay documented replication candidates.

## Consequences and limits

- The released axes are **not** the frozen 18,164 axis and do not replace it: identity,
  corpora and views keep the frozen axis, and `G_2000`/`G_5000` are ordered index subsets.
- The axis is a **selection on measurements**, not a biological claim. It says which genes
  the model should read and predict, and `rank_sum` deliberately erases effect magnitude.
- Seven sources make `pN ∈ 0..7`; the statistic is coarse by construction, and a gene that
  appears in two categories already carries most of the achievable evidence. The histogram
  and the per-category composition table are reported so the axis is not read as if
  `pN = 7` were a continuous quality score.
- Drug classes pool heterogeneous drugs inside one mechanism; the class-level DEG buys
  mechanism coverage at the price of washing out drug-specific effects.
- The `GSE270788` and Reichart controls are shared between categories, so category ranks are
  not statistically independent. This is declared, not corrected, and downstream comparison
  must not treat the seven sources as seven independent experiments.
- Cell-type labels for several candidates live in large serialized R objects; a source whose
  labels cannot be extracted at acceptable cost stays a gap rather than being re-clustered
  ad hoc.

## See

- [External data provenance and the Stage 2 step layout](0030-external-data-provenance-and-step-layout.md)
- [Stage 2 accepted products and independently reproducible figures](0062-stage2-accepted-products-and-reproducible-figures.md) — amended here for this step
- [Shared Stage 2 bases and indexed gene identity](0061-stage2-shared-bases-and-indexed-corpora.md) — corpus semantics amended here
- [Prepare the Stage 2 gene-axis decision](0065-stage2-gene-axis-decision-preparation.md) — superseded by this record
- `docs/research/stage2-gene-axis-cm-source-survey.md` — the source survey behind the roster
- `stage2/01_bases/AGENTS.md`, `stage2/01_bases/README.md`
- `GLOSSARY.md` — the gene-axis, corpus and view vocabulary
