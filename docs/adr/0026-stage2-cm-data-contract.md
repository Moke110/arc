# Stage 2 CM data — frozen corpus and filtered views

**Status:** accepted (2026-10-20); amended by [ADR 0033](0033-stage2-dataset-corpus.md) (2026-11-04)

Stage 2 needs cardiomyocyte perturbation data for model adaptation, but Stage 1's
`Union finetune data` is unusable for it: 978 landmarks, per-plate quantile
alignment to the L1000 pretrain, HG/`drug_idx` filtering, and a single collapsed
`cardiomyocyte` cell identity are all XPert-specific. `stage2/02_cm_data/` therefore
builds a model-agnostic corpus from the same four sources (`dtoxs_primary`,
`dtoxs_ipsc`, `screen_ctox` = GSE244740, `temp_env` = GSE262419) in **six
decisions**:

1. **Gene space = four-source intersection, latest HGNC symbol.** The four
   sources use symbols (38,478), Ensembl IDs (44,154), and `SYMBOL_entrez` rows
   (22,537 → 19,703 symbols), so a shared space requires an identity rule.
   Genes are matched alias/previous-symbol-aware and Entrez-verified where
   available; unresolvable or ambiguous genes are dropped and reported. The
   intersection is capped by the TempO-Seq panel (~18.7k), so `temp_env` gates
   the universe. Chosen over per-file native spaces because a frozen shared
   space is what makes cross-view and cross-dataset comparison reproducible
   without re-deciding gene identity at every adapter.
2. **`X` = raw counts**, duplicate probes summed by symbol before any
   transform; one `log1p_cp10k` convenience layer. No quantile alignment, batch
   correction, or HVG selection — those are model-specific and would make the
   pharmacology-tier comparison depend on a preprocessing artifact.
3. **Controls are rows**, with `control_group = dataset::plate` and a `vehicle`
   column, plus per-group means in `uns`. CPA-style models need a control
   population; mean-based models need the mean; one representation cannot serve
   both without discarding information.
4. **Drug identity = desalted parent InChIKey, stereochemistry retained**
   (`drug_key`). Connectivity-only keys are unusable: doxorubicin and epirubicin
   are 4′-epimers and share the same InChIKey connectivity block. Salt and
   free-base entries stay separate rows but share identity, so `Tamoxifen` /
   `Tamoxifen citrate` and `Estradiol` / `17beta-Estradiol` cannot leak across a
   cold-drug split. Structure wins over name: the DToxS `AMI` arm carries a
   pentyl-benzofuran homolog (CID 16218975), not amiodarone, so it enters
   `all-smiles` but not `atc-moa`.
5. **Three pharmacology tiers for `temp_env`**, built as filtered projections of
   one corpus: `atc-moa` (113 drug keys / 1,047 samples), `atc-moa-mesh`
   (252 / 2,310), `all-smiles` (460 / 4,206). MeSH is a separate tier because it
   admits herbicides, insecticides, carcinogens and disinfectants — chemical
   taxonomy, not pharmacology. The two curated drug panels use `atc-moa` with a
   cited **manual supplement** for real drugs missing from ChEMBL's mechanism
   table (currently SB-431542); the `temp_env` tiers stay strictly automatic,
   because that boundary *is* the comparison.
6. **The corpus is frozen without matrices.** `corpus/v1/` stores the gene
   set, the drug set, the per-source resolved gene sets (the drift baseline)
   and the exclusion report; no expression matrix is stored. Bundles re-read
   the sources and apply the frozen universes.
   The build **hard-fails** on any gene or drug outside the corpus, so a given
   corpus version always yields the same rows. Promoting a corpus version is a
   deliberate edit. Loaders are copied from
   `stage1/xpert-fda-miri-v3/03_cm_finetune_data/_datasets.py` (commit `c921f7d`)
   rather than imported, because `stage2/AGENTS.md` treats Stage 1 as
   non-dependency and Stage 2 must diverge on naming, dose provenance, and gene
   identity.

**Naming.** The two GEO sources use a short Stage 2 id: `screen_ctox` for
GSE244740 (ScreenSeq hiPSC-CM drug cardiotoxicity panel) and `temp_env` for
GSE262419 (TempO-Seq hiPSC-CM environmental chemical panel). The accession
remains the provenance anchor: it stays on the raw folder and on the read-only
Stage 1 baseline, is defined once in `cmdata/paths.py` (`SOURCE_ACCESSIONS`),
and is repeated in `stage2/02_cm_data/README.md`, `corpus/v1/README.md` and the
bundle `qa_report.md`.

**Update (2026-10-21) — physical plate ids.** GSE244740 numbers its
Optimization and Screening plates independently (`Optimization plate 1` and
`Screening plate 1` both appear as `plate 1`), but the series matrix titles were
being parsed without the assay token, so the loader collapsed each pair into one
`plate<n>`. `control_group = dataset::plate` therefore mixed two physical plates
for `screen_ctox::plate1..3`, and the per-group `control_mean` was a blend of the
two batches: the Screening baseline shifted by RMS 0.05–0.07 and 42-panel
contrasts were inflated (e.g. Amiodarone 470 → 361 significant genes with the
correct same-plate controls). Plate ids now carry the assay
(`optimization_plate1` … `screening_plate3`), so one `control_group` is one
physical plate. This is a bundle rebuild only: the corpus gene/drug universes
are unchanged, and no downstream training had consumed `datasets-1` yet.

**Consequences.** `datasets-1/` ships five views: `dtoxs_cm.h5ad` (957 + 118
rows, 46 drugs), `screen_ctox.h5ad` (2,520 + 552, 42), and three `temp_env` views
(1,047/2,310/4,206 perturbed rows, 113/252/460 drugs, each + 272 controls). The
intersection drops ~12k unmapped `dtoxs` identifiers and ~1.4k TempO-Seq genes
that are not shared; a future model needing them must rebuild with a new corpus
version, not silently widen a view. Stage 1 artifacts remain untouched and serve only as
the well-set regression baseline for QA.

**Update (2026-11-04) — corpus layer superseded in part.** ADR 0033 folds this
step into `stage2/02_dataset_corpus/`. What changes here: the corpus now carries
matrices, on each source's **native** gene axis, so decision 1's four-source
intersection survives as the **CM gene axis** (`genes.csv`, the projection target
for a benchmark or CM-model view) rather than as a corpus universe; decision 6
("the corpus is frozen without matrices") is retired; `datasets-1/` and its five
views are no longer produced by a data step — they belong to
`03_models_cmbench/` (frozen, to be rebuilt when the benchmark is redone), with
the 2026-10 build preserved read-only as `stage2/_frozen/cm-data-v1/`; and the
short ids `screen_ctox`/`temp_env` become `sctox`/`tenv`, with `dtoxs_cm`
splitting back into `dtoxs_primary` + `dtoxs_ipsc`. What stands: raw-count `X`
with one `log1p_cp10k` layer (decision 2), controls as rows with
`control_group = dataset::plate` and vehicle (3), the desalted-parent-InChIKey
drug identity and the structure-beats-name rule (4), the three pharmacology
tiers as rules rather than pre-applied filters (5), and the manual-supplement
boundary on the curated panels. The v1 gene/drug universes (18,164 symbols,
525 drug keys) are unchanged.
