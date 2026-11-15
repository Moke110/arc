# Stage 2 shared foundations

**Accepted (2026-11-13, ADR 0066).** This step's `output/` tree is the post-rebuild
foundation: `raw_counts/` -> `corpus/` -> `view-2k/`/`view-5k/`, with
`disease-deg/` and the released `gene-axis/` (`candidate_rank_sum.csv`,
`gene_axis_2000.csv`, `gene_axis_5000.csv`, `axis_composition.csv`,
`source_topn_coverage.csv`, `gene_axis_fill_rule.docx`, `manifest.json`). The measured truth, the successor envelope, the diagnostic PR
figure and the superseded build chain are **deleted**; git retains the code. Per
ADR 0066, `log/` is execution-run evidence rather than a change history: it was
deleted and regenerated at acceptance. [`log/acceptance.md`](log/acceptance.md)
is the acceptance record — replay results, digests, resource peaks, the "no
retired path read" evidence and the old-to-new mapping. The per-slice sections
below keep each ticket's contract, now as accepted science.

## Rebuild foundation: `raw_counts/` and `corpus/` (ADR 0066, #186)

The rebuild's first slice has landed. These two product areas are the new
`source/` -> axis chain; the superseded pre-0066 products these areas replaced were deleted by the
acceptance/removal ticket (#192); see *Removals* below.

| path | contents | reader |
|---|---|---|
| `output/raw_counts/` | `tenv.h5ad` (4,510 x 19,613 native int32), `tenv_probe_mapping.csv`, `tenv_axis_mapping.csv` (frozen `gene_index` -> native column), `tahoe.h5ad` (67,018 x 18,164 int32 CSR), `token_mapping.csv` (original Ensembl/token per frozen gene), `manifest.json` | `bases.raw_counts.open_tenv` / `open_tahoe` |
| `output/corpus/` | `tenv.h5ad` (4,478 rows = 4,206 treated + 272 vehicle), `tahoe_bulk.h5ad` (66,718 rows = 65,318 treated + 1,400 `DMSO_TF`), `exclusions.csv`, `manifest.json` | `bases.corpus.open_corpus` |

Raw counts are integer only: no normalization, no delta and no axis subsetting.
The corpus is the same counts on the frozen 18,164-gene axis with **controls as
ordinary rows**. Every row carries `role`, the source profile/well identity,
plate, cell line, desalted drug-parent InChIKey, `dose_um` and the original dose,
the source item and its item-manifest digest, and `matched_control_group`.
Matched controls are the same-plate DMSO wells for `tenv` and the same-plate,
same-cell-line `DMSO_TF` rows for `tahoe_bulk`; a treated row without a compatible
control is excluded and the reason is recorded in `exclusions.csv`. No
normalization, no delta, no filter and no HVG subsetting happens here — that is a
view's job.

```python
import sys; sys.path.insert(0, 'stage2/01_bases/output/api')
from bases.corpus import open_corpus
with open_corpus('stage2/01_bases/output/corpus') as corpus:
    tenv = corpus.source('tenv')
    counts = tenv.counts(0, 64)                # dense int32 window; never the whole matrix
    info = tenv.obs.iloc[42]                  # role, plate, dose, matched_control_group, ...
    controls = tenv.matched_control_rows(42)  # vehicle row positions for that group
```

Readers are the only supported access path. `CorpusSource.counts(start, stop)`
returns a dense int32 window and `batches(batch_size)` iterates it, so the full
observation x 18,164 matrix is never materialized. `control_groups()`,
`control_rows(group)` and `matched_control_rows(index)` resolve the control
correspondence the view builder consumes. `RawCounts` exposes the native/projected
payload, the probe/token evidence and the bridge frames.

### Build, replay and provenance

```bash
python stage2/01_bases/execution/build_raw_counts.py publish   # candidate -> verify -> replace
python stage2/01_bases/execution/build_raw_counts.py verify    # rebuild in scratch, compare every digest
python stage2/01_bases/execution/build_corpus.py publish
python stage2/01_bases/execution/build_corpus.py verify
```

Per ADR 0062 every build writes a candidate under `execution/checkpoints/` (never
inside `output/`), verifies it and only then replaces the accepted path; a failed
build leaves the accepted release usable. `verify` replays the build from the
pinned inputs and compares the manifest and every output digest, so building
twice produces identical digests. The measured acceptance runs were raw-count
publish/verify 16.8 s / 16.5 s and corpus publish/verify 1 m 9 s / 1 m 11 s, each
under 0.5 GiB peak RSS.

`raw_counts/manifest.json` carries the recipe, the frozen axis digest, every input
path + sha256 (no absolute paths), every output digest and the population counts.
Its provenance is explicit: the default publication mode **rehomes** the already
source-accepted tenv/Tahoe payloads, re-checks their bytes against the tracked
[`log/raw_counts_provenance.json`](log/raw_counts_provenance.json) record (the two
pre-0066 source-replay acceptance records and the native tenv build manifest,
restored verbatim from git history and digest-pinned by the manifest), re-derives
the Tahoe token mapping from the immutable source `gene_metadata.parquet`, and
records `provenance.mode: source-accepted-replay` with a pointer to
`acceptance.<tenv|tahoe>` in that record. Every manifest `inputs` key now resolves
in a merged checkout: the two source item manifests, the source gene metadata,
the frozen axis/crosswalk and the tracked provenance record — never the retired
`log/records/` or `output/{tenv,tahoe}-raw-counts/` paths. The stronger mode
`build_raw_counts.py --full-rebuild publish|verify` rebuilds both payloads from the
immutable `source/` trees through the original builders and replays them from
source (recorded cost: Tahoe build 1 h 35 m, independent source replay 2 h 35 m
under the 10 GiB cap); it stays available and is the mode to use if the rehomed
bytes are ever in doubt.

The compact manifests, the exclusion table and the mapping evidence
(`tenv_probe_mapping.csv`, `tenv_axis_mapping.csv`, `token_mapping.csv`) are
tracked; the multi-hundred-MB count payloads are gitignored and registered in
`docs/data-inventory.md`.

### DEG-table/rank contract (for tickets #187-#190)

`output/api/bases/deg_rank.py` publishes once the pure functions the disease
tickets import: `SymbolResolver`/`align_table` (symbol -> frozen `gene_index`;
ambiguous, off-axis and unmapped symbols are recorded, never guessed),
`deg_effect` (`|log2FC| * -log10(padj)`; `NA` padj is non-significant),
`per_gene_deg` (multi-contrast maximum collapse) and `dataset_deg_rank`
(within-dataset rank, unmeasured = `n_measured + 1`,
`rank_norm = rank / n_measured`, `gene_index` tie-break). `bases.gene_axis`
keeps the `deg_rank` re-export contract the accepted gene-axis evidence and
`test_gene_axis.py` assert by identity.

### Disease DEGs: `nHF`, `ICM` (ADR 0066, #187; `AF` removed by ADR 0069)

The disease side of the gene axis has landed for the two ticket-187 human
cardiomyocyte-derived categories. The `AF` category was removed by ADR 0069
because its only dataset (`GSE255612`) carried no ordering evidence (all
measured effect scores tied); the item and its evidence still exist outside 01,
but 01 no longer declares or publishes it. One product directory per category under
`output/disease-deg/` holds the compact Control-vs-Disease table and the
within-dataset rank of every source dataset (the level-1 inputs the fusion
ticket #190 consumes), plus the category catalog and release manifest:

| path | contents | reader |
|---|---|---|
| `output/disease-deg/nhf/` | `gse270788__nicm.*`, `reichart-cm__dcm_genotype_negative.*`, `catalog.json`, `manifest.json` | `bases.disease_deg.DiseaseDeg` |
| `output/disease-deg/icm/` | `gse270788__icm.*`, `gse270788__ami.*`, `catalog.json`, `manifest.json` | `bases.disease_deg.DiseaseDeg` |

`<dataset>.deg.csv` is one row per measured axis gene (sorted by the
within-dataset rank) with `gene_index`, `gene_symbol`, `contrast`, `log2fc`,
`stat`, `p_value`, `padj`, `base_mean`, `direction`, `significant`,
`n_contrasts`, `n_contrasts_supporting`, `effect`, `rank` and `rank_norm`.
`<dataset>.rank.csv` is the full frozen-axis frame `{gene_index, rank,
rank_norm}` (18,164 rows, unmeasured genes at `n_measured + 1` with
`rank_norm = (n_measured + 1) / n_measured`, per ADR 0066 *Rank
normalization* — never clamped to exactly 1.0). The effect score, maximum collapse, `NA`-padj and
tie-break rules are ticket #186's `bases.deg_rank` contract; the published
tables are written at 17 significant digits so their bytes encode the exact
order the rank was derived from.

```python
import sys; sys.path.insert(0, 'stage2/01_bases/output/api')
from bases.disease_deg import DiseaseDeg
db = DiseaseDeg('stage2/01_bases/output/disease-deg')
db.categories()                       # ['icm', 'nhf']
db.datasets('nHF')                    # ['gse270788__nicm', 'reichart-cm__dcm_genotype_negative']
table = db.table('nHF', 'gse270788__nicm')   # compact DEG table, rank-sorted
rank = db.rank('nHF', 'gse270788__nicm')     # 18,164-row frozen-axis rank frame
db.eligibility('ICM', 'gse270788__ami')      # tier/evidence/shared-control/time-window record
```

Every read is digest-checked against the category `manifest.json`. Each catalog
record pins its tracked inputs as repo-relative `inputs` and its item-local
derived artefacts as `derived_inputs` (`{relative, sha256}`, relative to the ADR
0030 item): the derived DEG/pseudobulk CSVs are gitignored per item, so their
paths are composed at build time from `source_item` + the relative name and are
never emitted as repo-relative literals (keeps `script/check_paths.sh` green in
a merged checkout). Every record also carries `n_case`/`n_control` for its
contrast (a bulk/pseudobulk DEG always has both arms). The four datasets and
their declared eligibility:

| category | dataset | evidence | tier | case/control | file |
|---|---|---|---|---|---|
| `nHF` | `gse270788__nicm` | CM pseudobulk, GSE270788 NICM arm vs shared Donors | 1 | 6 / 7 | `nhf/gse270788__nicm.*` |
| `nHF` | `reichart-cm__dcm_genotype_negative` | CM pseudobulk, Reichart genotype-negative DCM vs normal | 1 | 8 / 18 | `nhf/reichart-cm__dcm_genotype_negative.*` |
| `ICM` | `gse270788__icm` | CM pseudobulk, GSE270788 ICM arm vs shared Donors | 1 | 5 / 7 | `icm/gse270788__icm.*` |
| `ICM` | `gse270788__ami` | CM pseudobulk, GSE270788 AMI arm vs shared Donors (`time_window_basis: author_declared`) | 1 | 5 / 7 | `icm/gse270788__ami.*` |

Sources and derivations are ADR 0030 items: `disease_data/GSE270788` and
`disease_data/Reichart-CM` each carry `script/`
(download + offline verify + CM-gate/pseudobulk/DEG derivation), immutable
`source/`, declared `local/`, `author_log.json` and `manifest.sha256`. The
pseudobulk is summed per donor before testing (integer counts for the count
deposits, a per-donor mean CP10K for the normalised GSE255612 deposit), and the
one live DEG implementation is `output/api/bases/disease_deg.R` (limma: voom for
counts, `lmFit` for the log-scale deposit). The worker's `R_HOME` is resolved by
`output/api/bases/r_env.py`: an explicit argument or `ARC_R_HOME` wins, otherwise
the R of the conda env running the build (`<prefix>/lib/R`), otherwise the
`Rscript` on `PATH`. GSE270788's three arms share one
control set, so all three contrasts come from **one joint `~ 0 + arm` fit**
over the 23 donor pseudobulk profiles rather than three case-control fits.

**GSE270788 arm-size erratum (recorded in issue #187 and in the item author
log).** The ticket body and ADR 0066 say `Donor 6 / NICM 7`; the immutable
deposit says **Donor 7 / NICM 6 / ICM 5 / AMI 5** (both
`GSE270788_metadata.csv.gz` grouped by `sample`, and the series-matrix
`!Sample_characteristics_ch1` disease state). The deposit is authoritative: no
donor is relabelled, all 23 donors are used, and the seven Donor samples are the
shared control set. ADR 0066's roster counts are superseded by the deposit
(design unchanged): the exact evidence paragraph is in `log/acceptance.md`, and the
ADR erratum is to be landed as a separate docs-only ADR (0067) by the Lead — the
ADR text is deliberately not edited by this ticket.

```bash
python stage2/01_bases/execution/build_disease_deg.py publish --all   # candidate -> verify -> replace, per category
python stage2/01_bases/execution/build_disease_deg.py verify --all    # replay from pinned item DEGs, compare every digest
```

Both commands run under the arc python. The driver stages candidates in
`execution/checkpoints/`, verifies table/rank/catalog constraints and replaces
only after success; `verify` rebuilds and compares the manifest and every output
digest, so a replay is byte-identical. The compact CSVs, catalogs and manifests
are force-tracked (payloads otherwise gitignored); the item-local DEG tables and
pseudobulk matrices stay in the items' `local/` areas.

### Isogenic CM models: `genetic_cardiomyopathy` (ADR 0066, #188)

The hereditary-myocardial-disease category is the widest one: **eight datasets**
in one category at two tiers. The two Reichart genotype-split exports
(genotype-positive DCM donors and the 8 ARVC donors, tier 1, human CM
pseudobulk) and six isogenic hiPSC-CM bulk models (tier 2) that ADR 0066's
**Isogenic exception** keeps as their own evidence type beside the clinical
cohort. The catalog records the exception per dataset and at category level.

| path | datasets |
|---|---|
| `output/disease-deg/genetic_cardiomyopathy/` | `reichart-cm__dcm_genotype_positive.*`, `reichart-cm__arvc.*`, `gse224129__mybpc3_hcm.*`, `gse234141__ttntv_dcm.*`, `gse304575__lmna_h222p_dcm.*`, `gse305648__rbm20_dcm.*`, `gse253223__pkp2_knockdown_arvc.*`, `gse194297__dmd_cardiomyopathy.*`, `catalog.json`, `manifest.json` |

The compact tables and full-axis ranks follow the same contract as
`nHF`/`ICM`; what is new here is the **derivation** of the in-vitro bulk
sources and the way a multi-arm deposit becomes one dataset:

- **One shared derivation adapter**:
  `output/api/bases/bulk_deg.py` parses the GEO series matrix, binds matrix
  columns to samples by name (fail-closed bijection), names features, runs the
  one live DEG implementation (`disease_deg.R`) and writes the item-local
  matrix / sample sheet / long DEG table plus a tracked summary. Each item's
  `script/build_deg.py` only declares its design, so six items do not carry six
  copies of the source-adapter logic.
- **Integerisation is one shared, tested helper**
  (`bases.bulk_deg.integerise_counts`, round to nearest with ties to even, and a
  count-scale guard that refuses a normalised matrix declared as counts). It is
  used by `GSE224129` (RSEM expected counts) and `GSE253223` (estimated counts);
  the other four deposits are integer and take the strict-count voom path
  unchanged. The tracked summaries and the catalog `method` record the rule, the
  maximum absolute deviation, the library-size scale and the number of genes
  whose rounded row is all-zero (dropped by the count model, reported — never
  silently filtered). The shared `disease_deg.R` count contract is **not**
  relaxed: #186/#187/#189 recorded products depend on it.
- **Multi-arm deposits collapse by maximum effect or pool with a blocking
  covariate**, and the decision is recorded in both the item summary and the
  catalog: `GSE224129` (`+/-` and `-/-` vs `+/+`), `GSE234141` (four paired
  format × cell-type strata) and `GSE305648` (KO/PL/RQ vs WT, batch blocked)
  collapse per gene by the ticket-01 maximum-effect rule; `GSE253223` pools four
  harvest days in one fit with `day` as a blocking factor. Single-pair deposits
  (`GSE304575`, `GSE194297`) declare `pooled: false`.
- **The two Reichart datasets are derived read-only from ticket #187's tracked
  export** by `disease_data/Reichart-CM/script/build_genetic_deg.py`; the
  pseudobulk is not re-derived and #187's own tracked summary/table bytes are
  untouched. Catalog records state this in `method.derivation`.
- **Category catalog additions** on top of the `stage2-disease-deg-catalog-v1`
  shape: a category-level `isogenic_exception` record, a per-dataset
  `isogenic_exception` (`control_arm`, `why`) for the six models, and a
  per-dataset `strata` record (`decision`, `pooled`, `blocking_covariates`).
  Dataset records otherwise carry the mandatory
  `dataset_id, deg, rank, n_measured, n_significant, method, inputs` plus
  `derived_inputs`.

```bash
python stage2/01_bases/execution/build_disease_deg.py publish --category genetic_cardiomyopathy
python stage2/01_bases/execution/build_disease_deg.py verify  --category genetic_cardiomyopathy
```

The six isogenic items are built first (each item's own `script/fetch.sh` ->
`script/verify.sh` -> `script/build_deg.py`), then this driver publishes and
replays the category. The tests are
`execution/tests/test_disease_deg_isogenic.py` (integerisation contract, sample
binding, isogenic control selection, multi-arm pooling recorded and applied,
catalog shape, rank semantics, determinism, no absolute paths, genotype split).

### Drug-cardiotoxicity class sources (#189, ADR 0066)

`cardiotox_TOP2` and `cardiotox_other` treat a **mechanism class**, not a drug,
as evidence. Three already-owned items are read read-only: `disease_data/GSE243674`
(hiPSC-CM, DNR/DOX/EPI/MTX + trastuzumab), `pert_data/DToxS-iPSC`
(GSE217421 Conv arm; 54 codes) and `pert_data/GSE244740` (ScreenSeq; DOX/IDA/
Mitoxantrone).

| path under `output/disease-deg/` | contents |
|---|---|
| `cardiotox/` | the pinned crosswalk (`crosswalk.csv`), class selection (`classes.csv`), membership (`class_members.csv`), excluded codes (`excluded.csv`), per-source eligibility (`eligibility.csv`), rejected candidates (`rejected.csv`) and `manifest.json` |
| `cardiotox_TOP2/` | one pooled class DEG and within-dataset rank per source (`gse243674.*`, `dtoxs_ipsc.*`, `gse244740.*`), `catalog.json`, `manifest.json` |
| `cardiotox_other/` | one pooled class DEG and within-dataset rank per entered DToxS non-TOP2 class (`dtoxs_ipsc.atc_l4-L01EA/L01EB/L01EX.*`), `catalog.json`, `manifest.json` |

The compact tables (≈15.7 MB total) are force-tracked; the build's count and R
intermediates live only in `execution/checkpoints/` and are deleted after the
candidate is replaced. `bases.cardiotox.DiseaseDeg` and `.CategoryDeg` are the
readers; `catalog.json` carries `(category, dataset_id, class_id) -> deg, rank,
`n_measured`, `method`, `member_drugs`, `n_small_molecule_drugs` for the
fusion ticket, plus a top-level `catalog_sha256` = sha256 of the canonical
catalog JSON with that field removed (a file cannot contain its own digest); the
manifest repeats it. All three products use schema
`stage2-disease-deg-catalog-v1` with `datasets` as a list of full records.

**Pinned crosswalk (no hand-entered class).** Every `(source, code)` pair is
resolved to one annotation row in
`pert_data/Annotation/local/cm_compound_annotation.csv`, keyed on
`(lookup_name, dataset_inchikey)` where present: a source's own analogue (the
documented DToxS `AMI` pentyl-benzofuran homolog) keeps its own structure and is
never relabelled as the named drug. The class is the **ATC level-4
pharmacological subgroup** (five-character prefix); a code with several ATC
subgroups is assigned only when exactly one candidate is corroborated by its
`mesh_classes`, otherwise it is excluded as `ambiguous_multi_class`. A code the
annotation does not carry, or whose row has no ATC subgroup, is unresolved and
recorded. The partition is mandatory — one code, one class — so no well can enter
two class DEGs. Coverage of the three sources: 101 pairs (54 DToxS codes, 5
GSE243674 tokens, 42 GSE244740 compounds) → 69 resolved, 19 without an
annotation row, 2 without an ATC subgroup, 11 ambiguous; 36 classes carry a
resolved member.

**Class rule.** A class enters a source only with **at least three
small-molecule drugs** in that source (a small molecule is a structure-backed
annotation row; antibodies and peptides cannot satisfy it and are therefore
excluded — the ADR-declared trastuzumab, bevacizumab and rituximab are among the
unresolved codes). `cardiotox_TOP2` is the `L01DB` subgroup and requires every
member's `mesh_classes` to contain `Topoisomerase II Inhibitors`; a member
missing that evidence is excluded and recorded. Entered classes: the TOP2 class
in each of the three sources (4 / 4 / 3 drugs) and the DToxS `L01EA` (5),
`L01EB` (3) and `L01EX` (6) classes for `cardiotox_other`. The 76 excluded codes
and the classes below three are listed with their reason in `excluded.csv` /
`classes.csv`.

**Class DEG.** A class's drug wells are pooled into one Disease group and
compared with the same-plate/same-time vehicle in each source stratum, using the
source's own test: the DToxS deposit's pinned edgeR pipeline on its read counts
(TMM, trended dispersion, `exactTest`, `remove_small`, per cell line), DESeq2
`~ donor + group` on GSE243674's counts (3 h and 24 h), and the declared
two-sample Welch test on log2(CPM+1) for GSE244740 (one plate). The per-stratum
pooled contrasts then collapse per gene by the maximum `|log2FC| x -log10(padj)`
through `bases.deg_rank`, and the published table carries the winning contrast's
`log2fc`, `stat`, `p_value`, `padj`, `base_mean`, `direction`, `significant`,
`n_contrasts`, `n_contrasts_supporting`, `effect`, `rank`, `rank_norm`; the rank
file covers all 18,164 frozen genes with unmeasured genes at `n_measured + 1`.
Shared vehicle wells are declared per source in `eligibility.csv`; classes are
never treated as independent experiments. Each catalog record carries the real
pooled arm sizes: `n_case` is the pooled class-drug well count and `n_control`
the vehicle well count (the per-contrast arm size when it is constant, otherwise
the sum over the disjoint source strata), with `n_case_basis`/`n_control_basis`
recording the derivation and the per-stratum counts kept in `contrasts[]`."""


```bash
python stage2/01_bases/execution/build_disease_deg_cardiotox.py publish
python stage2/01_bases/execution/build_disease_deg_cardiotox.py verify
```

The build uses the repository's documented R build environment, the `biomni_e1`
conda env (R 4.5.3, edgeR 4.8.2, DESeq2 1.50.2; recorded as
`conda:biomni_e1`), because the `arc` env's R has neither package. Which
`Rscript` runs is configuration, not science: `output/api/bases/r_env.py`
resolves an explicit argument, then `ARC_RSCRIPT`, then the `Rscript` on `PATH`,
so exporting `ARC_RSCRIPT` to that env's `Rscript` is the documented setting on
this workstation. The pinned DToxS pipeline calls `decideTestsDGE`, removed in
edgeR >= 4; the R driver provides the documented two-group
`decideTests(method = "separate")` equivalent. The step-owned counts and R
outputs stay under `execution/checkpoints/`; the manifests record every input
path + sha256 with no absolute paths, and `verify` replays the build from the
pinned inputs and compares every digest.

This is the standalone foundation retained after the user removed the other
Stage 2 steps. Immutable raw inputs remain under `pert_data/`. Deleted historical
experiments and their checkpoints are outside this build. The layout follows
ADR 0062: reusable accepted products under `output/`, report figures and their
complete plotted values under `figure/`, reproducible execution/audit code and
tests under `execution/`, and every record under `log/`.

## Layout after the rebuild

```
output/
  api/            # the single public reader/metric implementation (ADR 0062)
  palette/        # shared figure palette
  fold_n_seed/    # fold/seed partitions (retained from the pre-0066 build)
  gene-identity/  # frozen 18,164 identity tables (genes.csv + tables/)
  raw_counts/     # source/ -> axis bridge: native tenv + projected Tahoe integers
  corpus/         # 18,164-axis integer counts, roles and matched-control groups
  tenv-hvg/       # retained pre-0066 descriptive and train-only rankings
  disease-deg/    # per-dataset DEG tables and within-dataset ranks (level-1 inputs)
  gene-axis/      # candidate_rank_sum.csv, gene_axis_{2000,5000}.csv,
                  #   axis_composition.csv, source_topn_coverage.csv,
                  #   gene_axis_fill_rule.docx, manifest.json
  view-2k/        # model-facing X/control/delta on G_2000
  view-5k/        # model-facing X/control/delta on G_5000
```

The retired pre-0066 layout (`shared-expression/`, `training-corpora/`,
`tenv-raw-counts/`, `tahoe-raw-counts/`, `tenv-truth/`, `tenv-truth-successor/`,
`historical-training-release/`) and the superseded `gene-axis/` file set are
**deleted** by acceptance ticket #192; git history retains the code. The old-to-new
mapping is recorded in [`log/acceptance.md`](log/acceptance.md).

## Identity

`output/gene-identity/` is the frozen 18,164-symbol axis every product shares:
`genes.csv` (the published `gene_index,symbol,hgnc_id,entrez_id,ensembl_id,match_type`
frame), `manifest.json` and the frozen `tables/` (`genes.csv` without `gene_index`,
`gene_crosswalk.csv`). `execution/build_identity.py` publishes and replays it;
`output/api/bases/identity.py` is the reader. `gene_index` is 0-based frozen row
order and is the repository's integer error-detection device (ADR 0061): an axis is
an ordered subset that keeps the original integers and is never renumbered.

## Retained tenv-HVG ranking

`output/tenv-hvg/` is the pre-0066 published product retained byte-identically:
`descriptive.csv` (the 4,206-treated-well per-well dispersion ranking that is the
tenv-HVG axis source) plus `fold_0/1/2.csv` (train-only **model-selection** rankings,
not an axis input) and `cohorts.csv`/`manifest.json`. `execution/build_hvg.py` and
`output/api/bases/hvg.py` remain its publisher/reader. Its pre-0066 input
`training-corpora/tenv.h5ad` is retired and is never read: `build_hvg.py`
`--training` now defaults to the live corpus payload `output/corpus/tenv.h5ad`
(published by `execution/build_corpus.py`) and refuses to run when that path is
absent. The builder still expects a training-style h5ad (`obs` with `sample_id`,
`drug_parent`, `dose_um`), so deriving that projection from the corpus is the
documented deferred re-point (ADR 0066 *Deferred*) and the accepted
`descriptive.csv` bytes and digests stay what the gene axis pins — the default no
longer names a retired payload.
## Released gene axes `G_2000` / `G_5000` (#190, ADR 0066)

`output/gene-axis/` is the **released** two-version gene axis: the master candidate
table, the two ordered axis files, their manifest, the per-category composition
table, the source top-`N` coverage table, the human-facing rule document and the
two-panel rank-sum evidence figure. The axes are built from the six
cardiomyocyte-derived disease-DEG categories (`output/disease-deg/`) fused at two
levels with the published tenv-HVG ranking. Nothing is recomputed here: the
disease-DEG ranks, the tenv-HVG ranking and the frozen 18,164-gene identity are
read-only inputs, and no corpus, view, model or truth is opened.

Public implementation: `output/api/bases/gene_axis.py` — the pure fusion/axis
definitions (`winning_tier`, `category_rank`, `fusion_rank_sum`,
`candidate_table`, `axis_fill`, `axis_intersection`, `axis_composition`,
`source_topn_coverage`, `source_topn_covered_counts`, `pn_distribution`,
`rank_sum_histogram`, `kneedle_knee`,
`second_difference_peak`, `detect_inflection`, `inflection_stability`) and the
digest-checked `GeneAxis` reader (`axis`, `candidate_table`, `composition`,
`coverage`). It also keeps the #186 re-export contract:
`SymbolResolver`, `complete_rank`, `deg_effect`, `per_gene_deg`, `rank_scores`,
`ALIGNMENT_CATEGORIES`, `PADJ_FLOOR`, `PADJ_THRESHOLD`, `normalized_rank` and
`AXIS_SIZE` are the same objects as in `bases.deg_rank`. Driver:
`execution/build_gene_axis.py publish|verify` (candidate-then-replace, replay
comparison of the manifest and every output digest).

**Two-level rank-sum fusion.** Within a category, the *winning tier* is the
category's minimum declared `tier`, and every dataset at that tier contributes
its within-dataset rank; ADR 0066's *Isogenic exception* additionally keeps a
record whose `isogenic_exception.applies` is true even above the winning tier,
so `genetic_cardiomyopathy` fuses its two tier-1 human CM-pseudobulk datasets
**plus** its six tier-2 isogenic hiPSC-CM models (winning tier 1, six records
retained by exception). Level 1 is the ascending `sum(rank_norm)` over those
datasets with `gene_index` tie-break; level 2 is
`rank_sum = sum(rank_norm)` over the **fused** category class ranks and the
tenv-HVG rank, **each normalized by the frozen-axis size** (`rank / 18164`;
ADR 0066 "each normalized by the frozen-axis size", GLOSSARY *Source rank*:
"tenv-HVG contributes its dispersion rank under the same normalization"). `pN`
counts how many of the fused sources place a gene in their top-`N` set;
`n_sources_measured` counts the fused sources that measure it.

**Non-informative sources (ADR 0068, AF removed by ADR 0069).** A dataset that
measures genes but whose continuous effect scores are **all tied** carries no
ordering evidence: its rank degenerates to the `gene_index` tie-break. The `AF`
category's only dataset, `gse255612__af`, was exactly that case — all 18,092
measured genes have `padj ≡ 1`, so every effect score was 0 — and ADR 0069
removed the whole category from 01 instead of re-publishing its table. The rule
remains implemented and tested: a non-informative dataset is not a fusion source,
and a category with no informative dataset is excluded from level-2 `rank_sum`,
`pN` and `n_sources_measured` (the manifest would list it under
`fusion.excluded_categories`). The released fusion declares five categories plus
tenv-HVG and every one is informative, so the manifest records
`fusion.n_sources: 6`, `fusion.excluded_categories: []` and
`release.pN_bounds: [0, 6]`.

**Axis fill rule.** `G_N` takes the top `0.5*N` tenv-HVG genes
(`tenv_hvg_rank <= N/2`; the rank is a complete `1..18,164` permutation, so this
is exactly `0.5*N` genes), then `0.5*N` more from the **remaining** candidates —
the step-1 genes are excluded, otherwise the 50/50 split breaks — by `pN`
descending, `rank_sum` ascending, `gene_index` ascending. The three keys mean:

- `pN` is the number of the **six fused sources** whose top-`N` set holds the gene
  (`N` matching the axis: `p2000` for `G_2000`, `p5000` for `G_5000`). It is
  **not** `n_sources_measured` (the number of sources that measure the gene):
  the two columns agree on only 2 of 18,164 rows at `p2000` and 96 at `p5000`, and
  substituting one for the other changes 289 / 425 axis genes.
- `rank_sum` is **smaller is better** (a descending value is better); sorting it
  the other way would change 300 / 1,004 axis genes.
- `gene_index` ascending is contractual because `rank_sum` is a sum of integer
  ranks and repeats often (among the remaining candidates at `N=2,000`, 2,038
  values are shared by 4,341 genes). For both released sizes the `0.5*N` boundary
  falls on a unique `(pN, rank_sum)` key, so the published axes do not depend on
  the tie-break; the published table is already in ascending-`gene_index` order,
  which makes a stable sort equivalent to the explicit key.

Both axes are ascending-`gene_index` subsets of the frozen space with the
original integers — never renumbered — and they are **not nested**:
`G_2000 ∩ G_5000 = 1,972` (Jaccard 0.3922; 98.60% of `G_2000` lies inside
`G_5000`, both axes carry genes the other does not). The manifest reports sizes,
intersection, per-category composition, input digests and the fill rule; no
consumer re-derives it.

Products (force-tracked compact CSVs; manifest sha256 `cfd3ba11…`):

| path | rows | digest |
|---|---:|---|
| `output/gene-axis/candidate_rank_sum.csv` | 18,164 | `ba328172…` |
| `output/gene-axis/gene_axis_2000.csv` | 2,000 | `cb07535f…` |
| `output/gene-axis/gene_axis_5000.csv` | 5,000 | `15e5c49e…` |
| `output/gene-axis/axis_composition.csv` | 40 | `6985ae92…` |
| `output/gene-axis/source_topn_coverage.csv` | 2 | `27e3a50d…` |
| `figure/data/gene_axis_rank_sum_histogram.csv` | 143 | `4ac96696…` |
| `figure/gene_axis_rank_sum_histogram.png` | — | `5c54608f…` |
| `figure/data/gene_axis_source_topn_coverage.csv` | 14 | `b405b5f1…` |
| `figure/gene_axis_source_topn_coverage.png` | — | `51ea4bc8…` |

The human-facing rule document `output/gene-axis/gene_axis_fill_rule.docx`
(sha256 `2bbe1f0e…`) opens with the source dataset information table (one row per
released dataset plus the tenv-HVG fusion source, with a footnote defining the
CM DEG count, the contrast semantics and the five non-obvious sample-size
bases), follows it with the audited-and-excluded
AF sources table (five rows: design/n, actual CM DEG count and one-clause
exclusion reason, closing on the declared-gap/re-admission statement), then
restates this fill rule with its five common misreadings and carries the source
coverage table and the evidence figure. `execution/make_gene_axis_docx.py`
regenerates it byte-identically from `source_topn_coverage.csv`, the PNG, the five
`output/disease-deg/<category>/catalog.json` files and
`execution/gene_axis_source_metadata.py` (the container is re-zipped with fixed
entry timestamps; `python-docx` is listed in `environment.yml`). It is a reading
aid, not a manifest output: no release consumer reads it.

The candidate table carries `symbol`, `gene_index`, `tenv_hvg_rank`, the five
`<category>_rank` columns, `rank_sum`, `n_sources_measured`, `p2000`, `p5000`
(all `pN` in `0..6`; `n_sources_measured` in `2..6`). `rank_sum` ranges 0.0237–5.8103 (median 3.0118).
Cross-source replication is coarse by design: over all candidates `p2000` is
10,125 genes at 0, 5,136 at 1, 2,046 at 2, 687 at 3, 141 at 4, 27 at 5 and 2 at
6; `p5000` is 3,942/5,298/4,355/2,808/1,298/402/61 for `pN` 0–6. Inside `G_2000`
the `pN` distribution is 563/580/687/141/27/2 for `pN` 1–6 (none at 0); inside
`G_5000` it is 653/782/1,804/1,298/402/61. The composition table publishes both
distributions so `pN` is not read as a continuous quality score.
Per-category composition of `G_2000` / `G_5000` (axis genes the source measures):
`nhf` 1,978 / 4,936, `icm` 1,908 / 4,751, `genetic_cardiomyopathy` 1,997 / 4,997,
`cardiotox_TOP2` 2,000 / 5,000, `cardiotox_other` 1,513 / 3,738; the tenv-HVG
halves are 1,000 / 2,500.

**Top-N source coverage.** Each fusion source's top-`N` set is the `N` genes
with `rank_s <= N` in the candidate table; the five active `<category>_rank`
columns and `tenv_hvg_rank` are complete `1..18,164` permutations, so each set
holds exactly `N` genes (for the five category sources every one of those genes
is measured; `tenv-HVG` ranks the full frozen axis), matching the `pN`
definition. Intersecting a top-`N` set with the released axis gives the coverage
below. The two sets carry the same number of genes, so a count is both the share
of the source's top-`N` the axis keeps and the share of the axis drawn from that
source's top-`N`; the sets overlap across sources (a gene can sit in several —
that is what `pN` counts), so the rows do not sum to 100%. For reference, a
random
`N`-gene set covers `N/18,164` = 11.01% (`G_2000`) / 27.53% (`G_5000`) of an
axis of the same size. The tables are descriptive, computed from
`candidate_rank_sum.csv` and `gene_axis_{2000,5000}.csv` only; they are not an
axis-selection criterion. The same numbers are the published
`output/gene-axis/source_topn_coverage.csv` (rows `N = 2000, 5000`; columns
`N, %tenv_hvg_gN, %nhf_gN, %icm_gN, %genetic_cardiomyopathy_gN,
%cardiotox_TOP2_gN, %cardiotox_other_gN`) and are reproduced in the rule
document.

`G_2000`, coverage of each source's top-2,000:

| source | top-2,000 | in `G_2000` | coverage |
|---|---:|---:|---:|
| `nhf` | 2,000 | 632 | 31.60% |
| `icm` | 2,000 | 651 | 32.55% |
| `genetic_cardiomyopathy` | 2,000 | 748 | 37.40% |
| `cardiotox_TOP2` | 2,000 | 578 | 28.90% |
| `cardiotox_other` | 2,000 | 712 | 35.60% |
| `tenv-HVG` | 2,000 | 1,174 | 58.70% |

`G_5000`, coverage of each source's top-5,000:

| source | top-5,000 | in `G_5000` | coverage |
|---|---:|---:|---:|
| `nhf` | 5,000 | 2,321 | 46.42% |
| `icm` | 5,000 | 2,350 | 47.00% |
| `genetic_cardiomyopathy` | 5,000 | 2,611 | 52.22% |
| `cardiotox_TOP2` | 5,000 | 2,224 | 44.48% |
| `cardiotox_other` | 5,000 | 2,416 | 48.32% |
| `tenv-HVG` | 5,000 | 3,275 | 65.50% |

`tenv-HVG` is the only source with guaranteed axis membership: the fill rule
places its top `N/2` (`1,000` / `2,500`) first and `174` / `775` more enter via
the `pN`/`rank_sum` half; the five category sources have no guaranteed gene.

**Source-coverage figure.** `figure/gene_axis_source_topn_coverage.png` is a
two-panel figure (2026-11-14 design request). Both panels draw the same layout:
one x group per fused source (`tenv-HVG` then the five categories; `af` is not a
source, ADR 0069) and two adjacent bars per group, one per released axis, the
bar's identity and colour encoding that axis (`G_2000` blue, `G_5000` orange);
groups keep a gap and the bars inside a group touch. Panel **a** ("Covered Gene
Count", y = covered genes) plots the number of that source's top-`N` genes the
axis contains — the numerator of the percentages above (`G_2000`
1,174/632/651/748/578/712; `G_5000` 3,275/2,321/2,350/2,611/2,224/2,416) — and
panel **b** ("Top-N Coverage Rate", y = coverage rate (%)) plots the same cells
as rates (`G_2000` 58.70/31.60/32.55/37.40/28.90/35.60; `G_5000`
65.50/46.42/47.00/52.22/44.48/48.32). The two values are
`|{g in G_N : rank_s(g) <= N}|` and `100 * count / N` from
`bases.gene_axis.source_topn_covered_counts` / `source_topn_coverage`; the
complete plotted values are `figure/data/gene_axis_source_topn_coverage.csv`
(both columns per bar, each panel record naming the column it draws), and the
CSV-only renderer `figure/script/render_gene_axis_source_topn_coverage.py`
reproduces the PNG byte for byte (sha256 `51ea4bc8…`) without opening `output/`.

**Evidence figure and inflection.** `figure/gene_axis_rank_sum_histogram.png`
is a two-panel figure, one panel per released axis (`N = 2000`, `N = 5000`;
2026-11-14 redesign). Both panels draw the identical gray candidate histogram
("gene candidates") over all 18,164 `rank_sum` candidates, and from every bar's
top downward the colored share of that bin the axis selects ("selected genes";
the colored total is `N`). No other plot-area element is drawn. Its complete
plotted values are `figure/data/gene_axis_rank_sum_histogram.csv`; the CSV-only
renderer `figure/script/render_gene_axis_rank_sum_histogram.py` reproduces the
PNG byte for byte (sha256 `5c54608f…`) without opening `output/`. The axes'
`rank_sum` intervals (`G_2000` 0.0237–4.9587, median 1.7507; `G_5000`
0.0237–5.0185, median 2.0439) stay in the manifest (`axes.<N>.rank_sum`) but are
no longer drawn. The Kneedle +
second-difference diagnostic is **reported as a diagnostic only**: Kneedle
(`S=1.0`, first interior arc maximum) puts a candidate at `N=1496`, the
second-difference peak at `N=18162`; they do not agree, so
`diagnostics.inflection = false` with its full parameters, and six
leave-one-source-out replays all reach the same decision (`stable: true`). It
never chooses `N`.

```bash
python stage2/01_bases/execution/build_gene_axis.py            # publish
python stage2/01_bases/execution/build_gene_axis.py --verify   # replay digests
python stage2/01_bases/figure/script/render_gene_axis_rank_sum_histogram.py \
    --csv stage2/01_bases/figure/data/gene_axis_rank_sum_histogram.csv \
    --png stage2/01_bases/figure/gene_axis_rank_sum_histogram.png
python stage2/01_bases/figure/script/render_gene_axis_source_topn_coverage.py \
    --csv stage2/01_bases/figure/data/gene_axis_source_topn_coverage.csv \
    --png stage2/01_bases/figure/gene_axis_source_topn_coverage.png
python stage2/01_bases/execution/make_gene_axis_docx.py        # rule document
```

`execution/tests/test_gene_axis.py` is fixture-only (8-gene fixture) apart from
four real-file replay groups (the rule document, its source-table/excluded-AF
contracts, the docx input-wiring checks, and the two figures): winning-tier
and isogenic-exception selection, level-1 single-dataset identity, level-1/2
sums and `gene_index` ties, `pN` boundaries 0 and 6 (and 0..7 when every declared
category is informative), degenerate-source exclusion from fusion/pN and the
`n_sources_measured`, `0.5*N + 0.5*N` fill,
pN/rank_sum/gene_index ordering, non-nested axes and intersection,
ascending-in-frozen-space axes, source top-`N` coverage (percentages, the empty
excluded column, out-of-contract ranks fail closed), the count table being the
percentages' exact numerators, histogram and inflection
stability, candidate-then-replace determinism (`--verify` replay), manifest
relative-path hygiene, the two-panel figure CSV-only replay and its
impossible-selection rejection, the two-panel source-coverage figure CSV-only
replay and its impossible-panel/value rejection, the rule document's
byte-identical replay, its source-table/excluded-AF contracts and its
sample-size footnote cross-check, the docx input wiring (`--catalog` and
`--metadata` are actually read), the `deg_rank`
re-export identity, and preservation of
the retired file set (55 tests). `execution/tests/test_gene_axis_source_metadata.py`
holds the metadata mapping's coverage, required-field, platform-evidence and
tenv-constant contracts (22 tests). Build/verify peak RSS ≈140 MiB under the 10 GiB cap. The pre-0066 `gene-axis/` file set (direction
composition, unweighted quota-free fused table, animal half-weight column,
tenv non-HVG curve, ortholog table, `inflection.json`) was retained on disk
byte-identically until #192's acceptance passed and is now deleted (see
*Removals*); nothing produced or consumed it.

## Model-facing views `view-2k/` / `view-5k/` (#191, ADR 0066)

`output/view-2k/` and `output/view-5k/` are the models' products: the integer
corpus normalized **once** on the full 18,164-gene library total (`CP10K`, then
natural `log1p`; ADR 0061's frozen normalization), then subset to one released
axis's genes in ascending original `gene_index`, with `control` = the mean of
that row's matched corpus controls normalized the same way and
`Δ = X − control` elementwise. The view **never** re-normalizes on its own gene
subset; there is one normalization and two gene subsets (ADR 0066 *Corpus and
view semantics*). No `{N}k_eval_truth/` directory is created here: the
evaluation redesign owns each view's truth later.

**Measured denominator basis (released axes, this build).** The tenv axis
carries **15.0548953%** (`G_2000`) and **37.4273748%** (`G_5000`) of the
full-axis CP10K mass, pooled over all 4,478 profiles (treated-only 15.0892505% /
37.3846128%); `tahoe_bulk` carries 10.8317768% / 24.7366390%. A view-local denominator would change the
per-gene input baseline by about 2.49× between the two axes and turn
`mean|Δ|` into an artifact of the denominator, which is why the one
normalization is on the **full axis, not the view**. ADR 0066's design-record
figures 12.9% / 36.3% came from an **interim pre-#190** axis; the measured value
recorded in each manifest is authoritative, and the decision is unaffected
(≈2.49× measured vs the ADR's implied ≈2.81×). Each manifest also records every
profile's full-axis library total through `sources.<source>.full_axis_library_total`
and `library_totals.csv`.

Layout (the step's matrix-and-metadata pattern, one h5ad per corpus source;
payloads gitignored, manifests force-tracked):

| path | rows × columns | bytes | sha256 |
|---|---|---:|---|
| `output/view-2k/manifest.json` | — | 5,754 | `44f0c0df…` |
| `output/view-2k/tenv.h5ad` | 4,478 × 2,000 | 59,912,368 | `c02f320c…` |
| `output/view-2k/tahoe_bulk.h5ad` | 66,718 × 2,000 | 1,312,787,850 | `8ed8655b…` |
| — rows | tenv 4,206 treated + 272 vehicle; tahoe_bulk 65,318 treated + 1,400 vehicle | | |
| `output/view-2k/library_totals.csv` | 71,196 | 3,727,979 | `dfaa5d7e…` |
| `output/view-5k/manifest.json` | — | 5,753 | `d7ac0227…` |
| `output/view-5k/tenv.h5ad` | 4,478 × 5,000 | 142,266,477 | `e1a18a9a…` |
| `output/view-5k/tahoe_bulk.h5ad` | 66,718 × 5,000 | 3,152,087,026 | `02abaed9…` |
| `output/view-5k/library_totals.csv` | 71,196 | 3,727,979 | `dfaa5d7e…` |

`du -sb`: `view-2k/` 1,376,433,951 B, `view-5k/` 3,298,087,235 B. Each h5ad
carries `X` (float32), `layers/control` and `layers/delta`; `obs` carries the
corpus row identity (`role`, plate/cell line, compound name, desalted parent
InChIKey/SMILES, source item + item-manifest digest, `dose_um`, source dose
column, `matched_control_group`) plus the row's full-axis `library_total`; `var`
carries the axis genes as original frozen integers (`gene_index`, `symbol`) in
ascending order, never renumbered to `0..N−1`. `manifest.json` records the
recipe, the full-axis normalization definition, input digests (corpus manifest +
both payloads + gene-axis manifest + axis file), output digests, per-source
counts, the measured basis and the environment; no absolute paths.

```python
import sys; sys.path.insert(0, 'stage2/01_bases/output/api')
from bases.view import open_view
with open_view('stage2/01_bases/output/view-2k') as view:
    tenv = view.source('tenv')
    x, control, delta = tenv.X(0, 64), tenv.control(0, 64), tenv.delta(0, 64)
    row = tenv.obs.iloc[42]                    # role, dose, structure, matched_control_group, library_total
    controls = tenv.matched_control_rows(42)   # vehicle row positions whose mean is `control`
    genes = tenv.gene_index                    # original frozen gene_index integers, ascending
```

Readers are the only supported access path. `View` / `ViewSource` expose
bounded dense windows (`X(start, stop)`, `control`, `delta`, `batches(batch_size)`
yielding `(start, X, control, delta)`), the carried `obs`/`var`, `gene_index`,
`library_totals()`, `control_groups()`, `control_rows(group)` and
`matched_control_rows(index)`, plus `digests()`; they never materialize a dense
observation × axis matrix and never open the corpus.

**Consumer contract for `02_models_tenv` (re-pointing deferred).** The model
step currently reads `output/training-corpora/tenv.h5ad` layers `control` and
`delta`; re-pointing it and re-running its accepted 2-fold × 2-seed matrix are
**deferred** (ADR 0066 *Deferred*) and were not attempted here. Its target
contract is `output/view-{N}k/tenv.h5ad`: `X` (log1p CP10K on the full-axis
denominator), `layers/control`, `layers/delta = X − control`, `obs` row identity,
`var.gene_index` ascending frozen integers. The view keeps **both corpus roles**
as rows (`role` is carried), so the model step selects `role == 'treated'`
(4,206 tenv rows; the retired training h5ad held treated rows only) instead of
treating vehicle rows as observations. The two axes share one scale, so the axis
comparison must still declare one hyperparameter protocol across axes or a
per-axis tuning budget (ADR 0066 *Tuning protocol for the axis comparison*).

```bash
python stage2/01_bases/execution/build_view.py publish   # candidate -> verify -> replace, both views
python stage2/01_bases/execution/build_view.py verify    # replay both from corpus + axis, compare every digest
```

Publication is candidate-then-replace per view; `verify` rebuilds both views
from the pinned corpus and gene-axis releases and compares the manifest and
every output digest (recorded reproducible byte-identically). Measured cost of
one publish run: 3 m 23.2 s, peak RSS 459,796 KiB; one verify replay:
3 m 32.5 s, peak RSS 440,348 KiB, both under the 10 GiB cap and streaming
256-row batches.
`execution/tests/test_view.py` is hermetic (17 tests): the **full-axis vs
view-local denominator** fixture asserting the full-axis result, `control` =
mean of the named matched controls, `Δ = X − control` exactly, replay
reproducing `Δ`, exactly 2,000/5,000 ascending frozen genes, no
`{N}k_eval_truth/`, determinism, tamper detection and no absolute paths. Run it
individually under the `arc` env.

## Readers

`output/api/bases/` is the single supported access path to every product. Import it
with `execution/public_api.py` (or `sys.path.insert(0, 'stage2/01_bases/output/api')`
from the repository root):

| product | reader |
|---|---|
| `raw_counts/` | `bases.raw_counts.open_tenv` / `open_tahoe` (bounded windows; probe/token evidence and the frozen-index bridge) |
| `corpus/` | `bases.corpus.open_corpus` (`counts`/`batches`, `matched_control_rows`, `control_groups`) |
| `view-2k/`, `view-5k/` | `bases.view.open_view` (`ViewSource` with `X`, `control`, `delta` windows and `library_totals`) |
| `disease-deg/` | `bases.disease_deg.DiseaseDeg` (`categories`, `datasets`, `table`, `rank`, `eligibility`; digest-checked per read) |
| `gene-axis/` | `bases.gene_axis.GeneAxis` (`axis`, `candidate_table`, `composition`, `coverage`, `digests`) |
| `gene-identity/` | `bases.identity` |

No reader ever materializes a dense observation x 18,164 matrix; every large payload
is streamed in bounded batches. Never open a payload h5ad with `backed="r"`: the
readers (or their bounded window/batch methods) are the only supported access path.

## Figure reproduction

`figure/gene_axis_rank_sum_histogram.png` is the released evidence figure (two
panels, one per axis; gray candidates plus the colored selected-gene overlay); its
complete plotted values are in `figure/data/gene_axis_rank_sum_histogram.csv` and
the CSV-only renderer is `figure/script/render_gene_axis_rank_sum_histogram.py`.
Re-rendering reads only that CSV (plus the pure `bases.palette` colour function)
and reproduces the committed PNG with the product payloads unreadable. The rule
document `output/gene-axis/gene_axis_fill_rule.docx` is generated by
`execution/make_gene_axis_docx.py` from the coverage CSV, the PNG, the five
`output/disease-deg/<category>/catalog.json` files and
`execution/gene_axis_source_metadata.py` (byte-identical replay; it reads
published products and one tracked metadata module, unlike the evidence
renderer). The retired
diagnostic PR figure (`figure/diagnostic_tenv_pr.*`) was deleted with the
measured-truth chain.

## Health checks and replay

```bash
script/check_paths.sh      # every external path literal in tracked files resolves
script/check_tracked.sh    # every executable file on disk is tracked
script/verify_all.sh       # every external item's offline manifest verification
```

`verify_all.sh` is the documented ~9.7 GiB page-cache job: run it alone, never beside
another heavy job. Product replay is per driver, all under the `arc` env:

```bash
python stage2/01_bases/execution/build_raw_counts.py verify       # ~15 s
python stage2/01_bases/execution/build_corpus.py verify           # ~1 m 11 s
python stage2/01_bases/execution/build_disease_deg.py verify --all
python stage2/01_bases/execution/build_disease_deg_cardiotox.py verify
python stage2/01_bases/execution/build_gene_axis.py --verify
python stage2/01_bases/execution/build_view.py verify             # ~3.5 m
```

Each `verify` rebuilds from the pinned inputs into `execution/checkpoints/`, compares
the manifest and every output digest, and never writes the accepted release. The
acceptance run's results, digests, peak RSS and the "no retired path read" evidence are
in [`log/acceptance.md`](log/acceptance.md).

Heavy work runs **one job at a time** through `execution/capped.sh` (systemd user
scope when available, otherwise an equivalent 10 GiB `RLIMIT_AS`) wrapped with
`/usr/bin/time -v`, and the measured peak is recorded. Disposable build state is
`execution/checkpoints/`, never inside `output/`.

## `log/`

ADR 0066 applies the repository rule strictly to this step: `log/` is execution-run
intermediate output and verification evidence, **not** a history of directory changes
(git is that history). It was deleted at the rebuild and regenerated at acceptance;
nothing is archived from the pre-0066 tree. [`log/acceptance.md`](log/acceptance.md)
and [`log/inputs.json`](log/inputs.json) are the step's acceptance record,
[`log/raw_counts_provenance.json`](log/raw_counts_provenance.json) carries the restored
pre-0066 source-replay evidence, and the per-slice replay records next to them
(`*_verify.json`) are the new build's evidence.

## Removals and the old-to-new mapping

Deleted after acceptance passed (ADR 0062/0066). Git retains the code.

| removed | replaced by |
|---|---|
| `output/shared-expression/` | `output/corpus/` (integer counts) + the view normalization |
| `output/training-corpora/` | `output/corpus/` + `output/view-2k/`/`view-5k/` |
| `output/tenv-raw-counts/`, `output/tahoe-raw-counts/` | `output/raw_counts/` (same bytes, one manifest) |
| `output/tenv-truth/`, `output/tenv-truth-successor/`, `output/historical-training-release/` | none; the evaluation redesign is deferred and each view will own its `{N}k_eval_truth/` |
| old `output/gene-axis/` file set | `candidate_rank_sum.csv`, `gene_axis_{2000,5000}.csv`, `axis_composition.csv` |
| truth/evaluation chain (`build_significance.py`, `bases.significance`, `bases.truth_build`, `bases.truth_successor`, `bases/evaluation`, `bases/metrics`, `bases/plots`, `deseq2_contrasts.R` + sealed copy, `evaluate_predictions.py`, `replay_pr.py`, `verify_hvg_metrics.py`, the diagnostic export/publication scripts, the diagnostic PR figure) | none; reworked with the deferred evaluation design |
| one-off migration/audit (`migrate_products.py`, `audit_metadata_bindings.py`, `verify_payload_digests.py`) | none; their job is done (records in git) |
| superseded build chain (`build_shared.py`, `build_corpora.py`, `verify_training_source.py`, `bases/corpora.py`, `build_native.py`, `prepare_tahoe.py`) | `build_raw_counts.py` (`publish`/`--full-rebuild`), `build_corpus.py`, `bases.raw_counts`, `bases.corpus` |
| `build_mouse_orthologs.R` + `mouse_human_orthologs.csv` | none; no animal DEG source survives eligibility |

## Recorded errata

**GSE270788 roster.** ADR 0066's category roster says "NICM arm (7 donors vs 6 Donor)"
and "the same six Donor controls". The immutable deposit disagrees: grouping
`GSE270788_metadata.csv.gz` by the `sample` column (the series-matrix
`!Sample_characteristics_ch1` agrees) measures **Donor 7 / NICM 6 / ICM 5 / AMI 5**, so
the shared control set is **seven Donor controls** and the NICM arm is 6 vs 7. The
deposit is used unchanged (the README's `nHF` table already records 6/7). The
`orig.ident` column is a single value `HF` across all arms and must not be used for
grouping. ADR text is not edited here; the erratum is to be landed as a docs-only ADR
(0067) by the Lead, with the evidence paragraph in `log/acceptance.md`.

**View basis.** ADR 0066 *Corpus and view semantics* quotes 12.9% / 36.3% of full-axis
CP10K for `G_2000`/`G_5000` from an interim pre-#190 axis. On the released axes
(after the #194 fusion amendment and both view rebuilds) the measured tenv shares are
**15.0548953% / 37.4273748%** (treated-only 15.0892505% / 37.3846128%); the 2.5x
denominator argument is unchanged and the manifests record the released-axis values as
authoritative. ADR 0066's figures are superseded, not re-based.

**Cosmetic, unfixed (out of scope).** `bases/disease_deg.py` emits a pre-existing
`SyntaxWarning: invalid escape sequence '\'` from a #187 docstring; it is cosmetic and
owned by a follow-up.

## Accepted results and limits

- Replay is byte-identical for every released product (see `log/acceptance.md`); the
  two axes are **not nested** (1,972 shared genes) and every axis index is in the
  frozen 18,164 space.
- `corpus/` carries no normalization, no delta and no HVG subsetting; `view-*` normalizes
  once on the **full-axis** library total and then subsets, so `X`/`control`/`Δ` share one
  scale across axes (ADR 0066).
- The axis is a selection on measurements, not a biological claim; `pN ∈ 0..6` is coarse,
  and the six fused sources are not statistically independent (shared `GSE270788` and
  Reichart controls, declared not corrected).
- `02_models_tenv` still reads the retired `training-corpora/tenv.h5ad`; re-pointing it at
  `view-2k/`/`view-5k/` and re-running its accepted matrix are **deferred** (ADR 0066
  *Deferred*), as is the evaluation redesign and recomputing the tenv-HVG ranking from the
  corpus (the default `--training` already names that live payload).
- Cell-type labels for some candidate sources live in large serialized R objects; those
  sources remain declared gaps rather than being re-clustered ad hoc.

## Historical boundary

The pre-0066 measured truth, successor envelope, diagnostic PR chain and the
unweighted quota-free/direction-composition axis evidence are retired science: do not
resurrect their semantics. Their ADRs (0061/0062/0063/0065) and the research notes under
`docs/research/` keep their original bytes as history; this README describes only the
post-rebuild tree.
