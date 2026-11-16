# Stage 2 shared foundations

**Accepted (2026-11-13, ADR 0066); disease side re-published 2026-11-15 (ADR 0071,
#208, with the v2 `rejected[]` rationale strings corrected and re-published by #212).**
This step's `output/` tree is the post-rebuild
foundation: `raw_counts/` -> `corpus/` -> `view-2k/`/`view-5k/`, with
`disease-deg/`, `disease-fusion/` and the released `gene-axis/`
(`candidate_rank_sum.csv`, `gene_axis_2000.csv`, `gene_axis_5000.csv`,
`axis_composition.csv`, `source_topn_coverage.csv`, `gene_axis_fill_rule.docx`,
`manifest.json`). The measured truth, the successor envelope, the diagnostic PR
figure and the superseded build chain are **deleted**; git retains the code. Per
ADR 0066, `log/` is execution-run evidence rather than a change history: it was
deleted and regenerated at acceptance. [`log/acceptance.md`](log/acceptance.md)
is the acceptance record — replay results, digests, resource peaks, the "no
retired path read" evidence and the old-to-new mapping; its pre-v2 sections are
sealed history and its §10–§12 are the v2 republish gates. The per-slice sections
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

### DEG-table/rank contract (for tickets #187-#190, v2 in #208)

`output/api/bases/deg_rank.py` publishes once the pure functions the disease
tickets import: `SymbolResolver`/`align_table` (symbol -> frozen `gene_index`;
ambiguous, off-axis and unmapped symbols are recorded, never guessed),
`deg_effect` (`|log2FC| * -log10(padj)`; `NA` padj is non-significant),
`per_gene_deg` (multi-contrast maximum collapse), the v2 **pass-any-contrast**
`deg` flag (`padj < 0.05` **and** `|log2FC| >= 2` in at least one declared
contrast) and `dataset_deg_rank` (DEG-only within-dataset rank: DEGs `1..n_DEG`,
every other frozen-axis gene `n_DEG + 1`, `rank_norm = rank / n_DEG`,
`gene_index` tie-break). The schema label moved `stage2-disease-deg-v1 ->
stage2-disease-deg-v2`. `bases.gene_axis`
keeps the `deg_rank` re-export contract the accepted gene-axis evidence and
`test_gene_axis.py` assert by identity.

### Disease DEGs: six categories (ADR 0066 #187-#189; v2 by ADR 0071, #208)

The disease side publishes one product directory per category — `nhf`, `icm`,
`af`, `genetic_cardiomyopathy`, `cardiotox_TOP2`, `cardiotox_other`, plus the
`cardiotox/` crosswalk — holding the compact Control-vs-Disease table and the
DEG-only within-dataset rank of every source dataset (the level-1 inputs the
fusion consumes), plus the category catalog and release manifest. The `AF`
category is published again: ADR 0069/0070 removed it, ADR 0071 re-admits it.
Roster in this release: **27 published datasets** (18 before #202-#206, +5 new
bases from #202/#203/#206, +4 AF from #204) across the six categories:
`nhf` 6 sources, `icm` 3, `af` 4, `genetic_cardiomyopathy` 8, `cardiotox_TOP2` 3,
`cardiotox_other` 3.

| path | contents | reader |
|---|---|---|
| `output/disease-deg/nhf/` | six `*.deg.csv` + `*.rank.csv` pairs, `catalog.json`, `manifest.json` | `bases.disease_deg.DiseaseDeg` |
| `output/disease-deg/icm/` | three pairs, `catalog.json`, `manifest.json` | `bases.disease_deg.DiseaseDeg` |
| `output/disease-deg/af/` | four pairs, `catalog.json`, `manifest.json` | `bases.disease_deg.DiseaseDeg` |
| `output/disease-deg/genetic_cardiomyopathy/` | eight pairs, `catalog.json`, `manifest.json` | `bases.disease_deg.DiseaseDeg` |
| `output/disease-deg/cardiotox_TOP2/`, `cardiotox_other/` | three class sources each, `catalog.json`, `manifest.json` | `bases.disease_deg.DiseaseDeg` |

`<dataset>.deg.csv` is one row per measured axis gene (sorted by the
within-dataset rank) with `gene_index`, `gene_symbol`, `contrast`, `log2fc`,
`stat`, `p_value`, `padj`, `base_mean`, `direction`, `significant`, **`deg`**
(the authoritative pass-any-contrast v2 flag), `n_contrasts`,
`n_contrasts_supporting`, **`n_contrasts_deg`**, `effect`, `rank` and
`rank_norm`. Note that the published `contrast` label is still the
maximum-effect (winning) contrast per gene — the flag, not the label, is what
decides DEG membership under v2. `<dataset>.rank.csv` is the full frozen-axis
frame `{gene_index, rank, rank_norm}` (18,164 rows, non-DEG genes at
`n_DEG + 1` with `rank_norm = (n_DEG + 1) / n_DEG` — never clamped to exactly
1.0, and a `n_DEG = 0` source ranks every gene at 1 and is excluded as
non-informative). The effect score, maximum collapse, `NA`-padj and
tie-break rules are ticket #186's `bases.deg_rank` contract; the published
tables are written at 17 significant digits so their bytes encode the exact
order the rank was derived from.

Each category `catalog.json` is the canonical `stage2-disease-deg-catalog-v1`
shape (the label is unchanged; the record inside is v2): `category`,
`directory`, `evidence_type`, `tier`, `n_deg` (axis-resolved DEG-union size),
`deg2000` (`DEG2000_d` size), `datasets[]` (per-source `tier`, `evidence_type`,
`n_deg`, `n_deg_in_deg2000`, `informative` + `informativeness`, `n_measured`,
`n_axis_unmeasured`, `n_significant`, `contrasts`, `inputs`, `derived_inputs`),
`rejected`, and — for every published category — the computed `tiering` record
(`criterion: axis_resolved_deg_union`, `limit`, `stop_reached`, `stopped_at_tier`,
`tiers_entered`, `tier_union`) plus `tier_sources` (the declared per-tier source
ladder). Both are derived from the declared per-source tiers by
`bases.disease_deg.tier_ledger`, never hand-written.

```python
import sys; sys.path.insert(0, 'stage2/01_bases/output/api')
from bases.disease_deg import DiseaseDeg
db = DiseaseDeg('stage2/01_bases/output/disease-deg')
db.categories()                       # ['af', 'cardiotox_TOP2', 'cardiotox_other',
                                      #  'genetic_cardiomyopathy', 'icm', 'nhf']
db.datasets('nHF')                    # ['gse270788__nicm', 'reichart-cm__dcm_genotype_negative',
                                      #  'gse35108__dcm_ipsc_cm', 'gse165303__dcm_lv_bulk',
                                      #  'gse116250__hf_lv_bulk', 'gse141910__hf_lv_bulk']
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
contrast (a bulk/pseudobulk DEG always has both arms). The four #187 datasets and
their declared eligibility (the other 23 published datasets, from #188/#189/#202-
#206, are listed per category in the catalogs and in the fusion ledger below):

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

Two optional design-spec fields extend the same worker without changing the
defaults (#202): `"block"` names a sample-sheet column and puts one consensus
intra-block correlation from `duplicateCorrelation()` into both `voom()` and
`lmFit()`, for sources where one subject contributes repeated profiles
(`disease_data/Kuppe-MI`'s donor x zone pseudobulk: the four infarct-heart IZ
profiles against the four normal-heart profiles, with the donor's other zone
profiles in the same joint `~ 0 + zone` fit); `"filter"` applies the
`edgeR::filterByExpr` rule (`min_count`/`min_total_count`/`large_n`/`min_prop`)
as a gene filter, written out in the worker because the pinned `arc` R env ships
limma but not edgeR. The filter only removes genes (the all-zero-row drop is
always applied) and reports the dropped count in the `DEG_SUMMARY` line; the
count contract itself (integer, non-negative counts) is unchanged.

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

### Per-disease fusion: level-1 rank and `DEG2000` (#205, ADR 0071)

The disease half of the v2 axis is DEG-only. Each accepted per-category product
is read back, each source's published compact table is reduced to its v2 DEG-only
rank (DEGs `1..n_DEG` by descending effect with the `gene_index` tie-break; every
other frozen-axis gene `n_DEG + 1`; `rank_norm = rank / n_DEG`), and the disease
fuses its informative sources at level 1:

```
disease_rank_sum = Σ_source rank_norm
disease rank     = ascending position of the sum, ties by gene_index
DEG union        = genes that are a DEG in at least one informative source
DEG2000_d        = DEG union ∩ (disease rank <= 2,000)   # all of it when smaller
```

A source with `n_DEG = 0` has no `rank_norm` (the ADR 0071 denominator is zero)
and a source whose measured effects are all tied carries no ordering evidence;
both are **excluded from the sum and recorded with their reason** (ADR 0068), so a
`NaN` can never silently drop a gene from a disease. Because a DEG's normalized
rank is `<= 1` while a non-DEG's is `(n_DEG + 1) / n_DEG > 1` in every source, the
DEG union occupies exactly the first disease ranks and `|DEG2000_d|` is always
`min(|union|, 2,000)`.

| path under `output/disease-fusion/` | contents | reader |
|---|---|---|
| `catalog.json` | schema `stage2-disease-fusion-v1`, `rank_contract: stage2-disease-deg-v2`, `rank_semantics: deg-only-v2`, `limit`, the frozen-axis pin, and per disease: `tier`, `evidence_type`, `declared_source_count`, `informative_source_count`, `n_deg` (DEG-union size), `deg2000` (`DEG2000_d` size), the per-source ledger (`n_deg`, `n_deg_in_deg2000`, `tier`, `flag_source`), `excluded_sources` with reason and `not_published` declared sources | `bases.disease_deg.DiseaseFusion` |
| `<directory>.deg2000.csv` | `gene_index, disease_rank, disease_rank_sum`, one row per `DEG2000_d` gene in rank order (a header-only file for an empty union) | `bases.disease_deg.DiseaseFusion.deg2000` |
| `manifest.json` | schema `stage2-disease-fusion-release-v1`, candidate-then-replace manifest (inputs, outputs, digests, `counts` = 6 diseases / 6 published / 7,088 `DEG2000` rows) | same |

The released ledger (per disease: sources / tiers entered / DEG union / `DEG2000`)
is `nHF` 6 / 1→2→3 / 1,661 / 1,661, `ICM` 3 / 1 / 1,647 / 1,647, `AF` 4 / 1→3 /
1 / 1 (`SLC9A7`), `genetic_cardiomyopathy` 8 / 1→2 / 1,482 / 1,482,
`cardiotox_TOP2` 3 / 2 / 4,514 / 2,000 and `cardiotox_other` 3 / 2 / 297 / 297.
Only `cardiotox_TOP2` reaches the 2,000 cut; the others publish their whole union,
which the build rule permits. `DEG2000` totals 7,088 rows over 5,282 distinct
genes. ADR 0071's drafted ledger carried 4,363 / 279 for the two cardiotox
classes; the published values are larger because the v2 flag is
pass-any-contrast rather than the collapsed winning-contrast row (see
[`log/acceptance.md`](log/acceptance.md) §10.3, which reproduces the old numbers by
dropping the `deg` column).

Public implementation: `output/api/bases/disease_deg.py` —
`source_deg_rank`, `level_one_rank` (also re-exported as
`bases.gene_axis.category_rank`, the single live level-1 implementation),
`fuse_disease`, `deg2000_rows`, `validate_catalog`/`build_catalog` (the
`tier`/`evidence_type`/`n_deg`/`deg2000`/`tiering` category-catalog shape and its
validation), `tier_ledger` (the tier-stop rule and the `tiering`/`tier_sources`
record) and the `DiseaseFusion` reader.

```bash
python stage2/01_bases/execution/build_disease_fusion.py publish
python stage2/01_bases/execution/build_disease_fusion.py verify
```

The driver is roster-neutral: it fuses whatever `output/disease-deg/` has
accepted and lists the rest of `execution/disease_deg_sources.py`'s declaration as
`not_published`. Every published v2 table carries the authoritative pass-any
`deg` flag, so the ledger records `flag_source = published_deg_flag`; a table
without the flag (a pre-#208 v1 table) is refused rather than re-derived.

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
  disease-fusion/ # per-disease level-1 rank_sum + DEG2000 cut (ADR 0071 level 1)
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
## Released gene axes `G_2000` / `G_5000` (#190 / #207, ADR 0066 + ADR 0071)

`output/gene-axis/` is the **released** two-version gene axis: the master candidate
table, the two ordered axis files, their manifest, the per-category composition
table, the source coverage table, the human-facing rule document and the two-panel
rank-sum evidence figure. The disease side is the **v2 rule of ADR 0071**:
DEG-only, tiered sources collected per disease and fused at two levels, with
`f2000` (not `pN`) selecting the axis's second half. The schema label stays
`stage2-gene-axis-release-v1` and the manifest declares
`release.algorithm = f2000-deg2000-v2`; a pre-v2 manifest is refused by
`build_gene_axis.py --verify` with a self-describing error. Nothing is recomputed
here: the published `DEG2000` tables and disease ranks, the tenv-HVG ranking and
the frozen 18,164-gene identity are read-only inputs, and no corpus, view, model or
truth is opened.

Public implementation: `output/api/bases/gene_axis.py` — the pure v2 fusion/axis
definitions (`candidate_columns`, `category_rank`, `deg2000_frame`,
`f2000_count`, `cross_disease_rank_sum`, `candidate_table`, `axis_fill`,
`axis_intersection`, `deg2000_composition`, `f2000_distribution`,
`source_topn_coverage`, `rank_sum_histogram`, `kneedle_knee`,
`second_difference_peak`, `detect_inflection`, `inflection_stability`) and the
digest-checked `GeneAxis` reader (`manifest`, `axis_sizes`, `digests`,
`candidate_table`, `axis_frame`, `axis`, `composition`, `coverage`). It also keeps
the #186 re-export contract: `SymbolResolver`, `complete_rank`, `deg_effect`,
`per_gene_deg`, `rank_scores`, `ALIGNMENT_CATEGORIES`, `PADJ_FLOOR`,
`PADJ_THRESHOLD`, `normalized_rank` and `AXIS_SIZE` are the same objects as in
`bases.deg_rank`, and the v1 names (`winning_tier`, `fusion_rank_sum`,
`pn_distribution`, `source_topn_covered_counts`, the `p2000`/`p5000`/
`n_sources_measured` columns) are **removed** — `test_gene_axis.py` asserts their
absence so they cannot creep back. Driver: `execution/build_gene_axis.py
publish|verify` (candidate-then-replace, replay comparison of the manifest and
every output digest).

**Two-level rank-sum fusion (ADR 0071).** Level 1 runs inside a disease: each
informative source contributes its DEG-only `rank_norm` (`1..n_DEG` by descending
effect, `gene_index` tie-break; every other frozen-axis gene `n_DEG + 1`;
`rank_norm = rank / n_DEG`), and the disease's rank is the ascending position of
`sum(rank_norm)`. Every declared source of the disease contributes — the tiers are
a *collection* order that decides when to stop adding sources, not a filter on the
fusion; `genetic_cardiomyopathy` therefore fuses its two tier-1 Reichart exports
**plus** its six tier-2 isogenic hiPSC-CM models, and `nHF` fuses its six sources
across tiers 1-3. Level 2 is `rank_sum = Σ_d rank_d / n_d` over the **six disease
terms only**, where `n_d` is the disease's axis-resolved DEG-union size and a gene
outside `DEG2000_d` takes `(n_d + 1) / n_d`. tenv-HVG is **not** a term (it
supplies the axis's first half), `rank_sum` is **ascending-is-better**, and
`f2000` = the number of diseases whose `DEG2000_d` contains the gene is the
disease-side selection signal that replaced `pN`.

**Non-informative sources (ADR 0068, kept by ADR 0071).** A dataset that measures
genes but whose continuous effect scores are **all tied**, or that has no DEG at
all (every gene then ranks at `n_DEG + 1 = 1`), carries no ordering evidence.
The rule remains implemented and tested: such a dataset is not a fusion source and
is recorded with its reason, and a disease with no informative source is excluded
from `rank_sum` and `f2000` (the manifest would list it under
`fusion.excluded_diseases`). In this release `AF` shows both faces of the rule:
three of its four sources (`gse255612__af`, `gse238242__af`, `gse269882__af`) have
0 v2 DEGs and are excluded with recorded reasons, while `gse138252__af` has one
(`SLC9A7`) and carries the whole category. The manifest records
`fusion.n_declared_diseases: 6`, `fusion.n_diseases: 6`,
`fusion.excluded_diseases: []` and `release.f2000_bounds: [0, 6]`.

**Axis fill rule.** `G_N` takes the top `0.5*N` tenv-HVG genes
(`tenv_hvg_rank <= N/2`; the rank is a complete `1..18,164` permutation, so this
is exactly `0.5*N` genes), then `0.5*N` more from the **remaining** candidates —
the step-1 genes are excluded, otherwise the 50/50 split breaks — that have
`f2000 >= 1`, ordered by `f2000` descending, `rank_sum` ascending, `gene_index`
ascending. The three keys mean:

- `f2000` is a count of **diseases** (`0..6`), not of sources and not of
  top-`N` sets: it counts how many diseases published the gene in their
  `DEG2000`. It replaced `pN` because `pN` counted a source's *top-`N` overlap*
  over a measured-gene rank — an object the DEG-only rule no longer produces.
- `rank_sum` is **ascending is better** (a lower value selects first): rank 1 is
  the strongest DEG and a gene outside `DEG2000_d` pays `(n_d + 1) / n_d > 1` in
  that term. Its absolute scale is dominated by the smallest-`n_d` disease —
  `AF`'s `n_d = 1` makes every non-`SLC9A7` gene pay exactly `2.0`, about `2×` the
  five large-`n_d` terms combined — which is why the evidence figure has one spike
  at ≈`7.0` (the `f2000 = 0` mass at the constant maximum). That is an artefact of
  AF's one-gene denominator, not a fusion failure: within an `f2000` stratum the
  AF term is a constant for every non-AF gene, so ordering is unaffected.
- `gene_index` ascending is the deterministic last key.

When fewer than `0.5*N` remaining candidates have `f2000 >= 1`, the shortfall is
filled from unused tenv-HVG genes in ascending tenv-HVG rank order — never from a
gene on disease evidence it does not have. The fallback is **implemented and
exercised by tests but 0 in this release** (4,906 / 4,330 eligible candidates for
1,000 / 2,500 slots; `axes.<N>.second_half.tenv_hvg_fallback = 0`).

Both axes are ascending-`gene_index` subsets of the frozen space with the
original integers — never renumbered — and they are **nested**: `G_2000 ⊂ G_5000`
(the manifest measures `|G_2000 ∩ G_5000| = 2,000`, Jaccard 0.4), because `f2000`
and `rank_sum` do not depend on `N`. The v1 `pN` fill was axis-matched
(`p2000` for one axis, `p5000` for the other), so its axes were *not* nested. The
manifest reports sizes, intersection, per-disease composition, input digests and
the fill rule; no consumer re-derives it.

Products (force-tracked compact CSVs; manifest sha256 `a344f9e2…`, 28,216 B):

| path | rows | digest |
|---|---:|---|
| `output/gene-axis/candidate_rank_sum.csv` | 18,164 | `4c1a98d9…` |
| `output/gene-axis/gene_axis_2000.csv` | 2,000 | `637bb8e1…` |
| `output/gene-axis/gene_axis_5000.csv` | 5,000 | `6f467c2c…` |
| `output/gene-axis/axis_composition.csv` | 32 | `51192501…` |
| `output/gene-axis/source_topn_coverage.csv` | 14 | `8db5fcbf…` |
| `figure/data/gene_axis_rank_sum_histogram.csv` | 143 | `c402554a…` |
| `figure/gene_axis_rank_sum_histogram.png` | — | `874f6020…` |
| `figure/data/gene_axis_source_topn_coverage.csv` | 16 | `00a04341…` |
| `figure/gene_axis_source_topn_coverage.png` | — | `37c24d17…` |

The human-facing rule document `output/gene-axis/gene_axis_fill_rule.docx`
(sha256 `e7fbc26a…`, 94,590 B) opens with the v2 build rule (`padj < 0.05` and
`|log2FC| >= 2`, pass-any-contrast; DEG-only source rank; level-1 fusion; the
tier-1→2→3 collection with the axis-resolved 2,000 stop; the `f2000`/`rank_sum`
fill and the tenv-HVG fallback), then the source dataset table (表 1: 28 rows =
27 published datasets + the tenv-HVG fusion source, with a footnote defining the
v2 CM DEG count, the contrast semantics and the ten non-obvious sample-size
bases), 表 2 (one row per disease: source-dataset count, tier path, DEG union,
`DEG2000` and a note stating the cut or the shortfall), the four misreading
clarifications, the coverage table in the published long format and the evidence
figure. `execution/make_gene_axis_docx.py` regenerates it byte-identically from
`source_topn_coverage.csv`, the PNG, the six
`output/disease-deg/<category>/catalog.json` files and
`execution/gene_axis_source_metadata.py` (the container is re-zipped with fixed
entry timestamps; `python-docx` is listed in `environment.yml`). It is a reading
aid, not a manifest output: no release consumer reads it.

The candidate table carries `symbol`, `gene_index`, `tenv_hvg_rank`, the six
`<category>_rank` columns (DEG2000 rank, blank outside `DEG2000_d`), `rank_sum`
and `f2000`; it carries no `p2000`/`p5000`/`n_sources_measured`. `rank_sum` ranges
2.6624–7.0054 and `f2000` over all candidates is 12,882 genes at 0, 3,835 at 1,
1,146 at 2, 248 at 3, 48 at 4, 5 at 5 and none at 6. Inside `G_2000` the `f2000`
distribution is 624/258/817/248/48/5/0 and inside `G_5000`
1,548/2,005/1,146/248/48/5/0 for `f2000` 0–6 (`f2000 = 0` genes are tenv-HVG-half
genes). `axis_composition.csv` publishes the declared halves, the (zero) fallback,
that distribution and the per-disease `DEG2000` coverage. Per-disease coverage of
`G_2000` / `G_5000` (`|DEG2000_d ∩ G_N| / |DEG2000_d|`): `nhf` 638/1,661 (38.41%)
and 1,151/1,661 (69.30%), `icm` 550/1,647 (33.39%) and 1,042/1,647 (63.27%),
`af` 0/1 and 1/1, `genetic_cardiomyopathy` 635/1,482 (42.85%) and 1,068/1,482
(72.06%), `cardiotox_TOP2` 831/2,000 (41.55%) and 1,740/2,000 (87.00%),
`cardiotox_other` 199/297 (67.00%) and 256/297 (86.20%).

**`DEG2000` coverage (the long-format table).** The v2 coverage table
`output/gene-axis/source_topn_coverage.csv` keeps its **v1 filename** for
compatibility with the pinned manifests — renaming it would re-pin the gene-axis
manifest, both view manifests and force a republish — but its content is the v2
`DEG2000` table, in **long format** with an explicit denominator per row:
`N, n_axis, group, group_label, numerator, denominator, rate, note` (14 rows: two
axes × seven groups). The numerator is `|DEG2000_d ∩ G_N|` for a disease and
`|top_N(tenv) ∩ G_N|` for `tenv_hvg`; the denominator is `|DEG2000_d|` for a
disease and `N` for `tenv_hvg`; `rate = 100 * numerator / denominator`. The v1
wide format (`%<source>_gN` columns) divided everything by `N`, which made AF's
one gene read as a 100% top-`N` coverage; the long format states each group's own
set size instead. Because the denominator is each group's own set, the rows are
shares of that set, not of the axis; a group with an empty `DEG2000` would keep an
empty row (ADR 0068) and none does in this release. The table is descriptive and
not an axis-selection criterion, and the same numbers are reproduced in the rule
document.

`G_2000`, coverage of each group's own set:

| group | numerator | denominator | rate |
|---|---:|---:|---:|
| `nhf` | 638 | 1,661 | 38.41% |
| `icm` | 550 | 1,647 | 33.39% |
| `af` | 0 | 1 | 0.00% |
| `genetic_cardiomyopathy` | 635 | 1,482 | 42.85% |
| `cardiotox_TOP2` | 831 | 2,000 | 41.55% |
| `cardiotox_other` | 199 | 297 | 67.00% |
| `tenv_hvg` | 1,079 | 2,000 | 53.95% |

`G_5000`, the same seven groups:

| group | numerator | denominator | rate |
|---|---:|---:|---:|
| `nhf` | 1,151 | 1,661 | 69.30% |
| `icm` | 1,042 | 1,647 | 63.27% |
| `af` | 1 | 1 | 100.00% |
| `genetic_cardiomyopathy` | 1,068 | 1,482 | 72.06% |
| `cardiotox_TOP2` | 1,740 | 2,000 | 87.00% |
| `cardiotox_other` | 256 | 297 | 86.20% |
| `tenv_hvg` | 3,068 | 5,000 | 61.36% |

`tenv-HVG` is the only group with guaranteed axis membership: the fill rule
places its top `N/2` (`1,000` / `2,500`) first and `79` / `568` more enter via the
`f2000`/`rank_sum` half; no disease has a guaranteed gene. `AF`'s `0/1` and `1/1`
are the honest single-gene denominator: its one `DEG2000` gene (`SLC9A7`) is in
`G_5000` but not `G_2000`, because every one of `G_2000`'s 1,000 disease-half slots
is an `f2000 >= 2` gene.

**Source-coverage figure.** `figure/gene_axis_source_topn_coverage.png` is a
two-panel figure (2026-11-14 design request; labels corrected in #210). Both
panels draw the same layout: one x group per fusion group (`tenv-HVG` then the six
diseases, `af` included) and two adjacent bars per group, one per released axis,
the bar's identity and colour encoding that axis (`G_2000` blue, `G_5000`
orange); groups keep a gap and the bars inside a group touch. Panel **a**
("Covered Gene Count", y = covered genes) plots the numerator (`G_2000`
1,079/638/550/0/635/831/199; `G_5000` 3,068/1,151/1,042/1/1,068/1,740/256) and
panel **b** ("Coverage of Each Group", y = share of group covered (%)) plots the
same cells over each group's own denominator (`G_2000`
53.95/38.41/33.39/0.00/42.85/41.55/67.00; `G_5000`
61.36/69.30/63.27/100.00/72.06/87.00/86.20). The v1 panel title was
`Top-N Coverage Rate`, which was untrue once the denominator became the group's own
set — AF's `1-of-1` cell would have read as a 100% "top-`N`" coverage. The values
are `|DEG2000_d ∩ G_N|` and `100 * numerator / denominator` from
`bases.gene_axis.source_topn_coverage`; the
complete plotted values are `figure/data/gene_axis_source_topn_coverage.csv`
(both columns per bar, each panel record naming the column it draws), and the
CSV-only renderer `figure/script/render_gene_axis_source_topn_coverage.py`
reproduces the PNG byte for byte (sha256 `37c24d17…`) without opening `output/`.

**Evidence figure and inflection.** `figure/gene_axis_rank_sum_histogram.png`
is a two-panel figure, one panel per released axis (`N = 2000`, `N = 5000`;
2026-11-14 redesign, x annotation corrected in #210). Both panels draw the
identical gray candidate histogram ("gene candidates") over all 18,164 `rank_sum`
candidates, and from every bar's top downward the colored share of that bin the
axis selects ("selected genes"; the colored total is `N`). No other plot-area
element is drawn. Its complete plotted values are
`figure/data/gene_axis_rank_sum_histogram.csv`; the CSV-only renderer
`figure/script/render_gene_axis_rank_sum_histogram.py` reproduces the PNG byte for
byte (sha256 `874f6020…`) without opening `output/`. The x annotation reads
`rank_sum (lower is better)`: ascending is better, and the pre-#210 label said the
opposite. The axes' `rank_sum` intervals stay in the
manifest (`axes.<N>.rank_sum`: `G_2000` 2.6625–7.0055, median 5.9610; `G_5000`
2.6625–7.0055, median 6.2573) but are no longer drawn. The single spike at
≈`7.0` is the `f2000 = 0` mass at the constant maximum of the six-term sum (AF's
`n_d = 1` term costs every non-`SLC9A7` gene exactly `2.0`); it is the rule's
arithmetic, not a rendering or fusion defect. The Kneedle +
second-difference diagnostic is **reported as a diagnostic only**: Kneedle
(`S=1.0`, first interior arc maximum) puts a candidate at `N=1061`, the
second-difference peak at `N=106`; they do not agree, so
`diagnostics.inflection = false` with its full parameters, and six
leave-one-disease-out replays all reach the same decision (`stable: true`; the
manifest's `stability.method` records `leave-one-source-out`). It never chooses
`N`.

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
four real-file replay groups (the rule document, its source-table/disease-table
contracts, the docx input-wiring checks, and the two figures): level-1
single-dataset identity, the `Σ rank_norm` sum with `gene_index` ties, the
non-DEG `(n_DEG + 1) / n_DEG` rule, zero-DEG and out-of-contract frames failing
closed, `f2000` counting diseases (not datasets) and `tenv` provably **not** a
`rank_sum` term, the union-size denominator, the `0.5*N + 0.5*N` fill with
`f2000`-desc-/`rank_sum`-asc ties and `gene_index` as the last key, the tenv-HVG
fallback (exact half boundary, unused-tenv order, exhaustion, and never admitting
a gene with no disease evidence), nested axes and the intersection, the six
declared categories and their release contract, `DEG2000` coverage carrying its
own denominator, the `f2000` distribution bound, histogram and inflection
stability, candidate-then-replace determinism (`--verify` replay), the pre-v2
manifest refusal, manifest relative-path hygiene, the two-panel figure CSV-only
replay and its impossible-selection rejection, the two-panel source-coverage
figure CSV-only replay and its impossible-panel/value rejection, the rule
document's byte-identical replay, its source-table and per-disease 表 2 contracts
and its sample-size footnote cross-check, the docx input wiring (`--catalog` and
`--metadata` are actually read), the retirement assertions (the v1 names are
absent and the retired gene-axis file set stays retired), and the `deg_rank`
re-export identity (67 tests).
`execution/tests/test_gene_axis_source_metadata.py`
holds the metadata mapping's coverage, required-field, platform-evidence and
tenv-constant contracts (25 tests). Build/verify peak RSS ≈140 MiB under the 10 GiB cap. The pre-0066 `gene-axis/` file set (direction
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
carries **13.1626914%** (`G_2000`) and **33.8379148%** (`G_5000`) of the
full-axis CP10K mass, pooled over all 4,478 profiles (treated-only 13.2072122% /
33.7895449%); `tahoe_bulk` carries 7.8913108% / 16.7625796%. A view-local denominator would change the
per-gene input baseline by about 2.57× between the two axes and turn
`mean|Δ|` into an artifact of the denominator, which is why the one
normalization is on the **full axis, not the view**. ADR 0066's design-record
figures 12.9% / 36.3% came from an **interim pre-#190** axis; the measured value
recorded in each manifest is authoritative, and the decision is unaffected
(≈2.57× measured vs the ADR's implied ≈2.81×). Each manifest also records every
profile's full-axis library total through `sources.<source>.full_axis_library_total`
and `library_totals.csv`.

Layout (the step's matrix-and-metadata pattern, one h5ad per corpus source;
payloads gitignored, manifests force-tracked):

| path | rows × columns | bytes | sha256 |
|---|---|---:|---|
| `output/view-2k/manifest.json` | — | 5,756 | `ac5a1a97…` |
| `output/view-2k/tenv.h5ad` | 4,478 × 2,000 | 58,627,248 | `ff821207…` |
| `output/view-2k/tahoe_bulk.h5ad` | 66,718 × 2,000 | 1,262,131,874 | `88739add…` |
| — rows | tenv 4,206 treated + 272 vehicle; tahoe_bulk 65,318 treated + 1,400 vehicle | | |
| `output/view-2k/library_totals.csv` | 71,196 | 3,727,979 | `dfaa5d7e…` |
| `output/view-5k/manifest.json` | — | 5,747 | `4620805d…` |
| `output/view-5k/tenv.h5ad` | 4,478 × 5,000 | 138,912,355 | `9c818600…` |
| `output/view-5k/tahoe_bulk.h5ad` | 66,718 × 5,000 | 3,012,415,752 | `5f16fdfb…` |
| `output/view-5k/library_totals.csv` | 71,196 | 3,727,979 | `dfaa5d7e…` |

`du -sb`: `view-2k/` 1,324,492,857 B, `view-5k/` 3,155,061,833 B; the payload
bytes are the v2 axes' (the axis gene set changed, so `X`/`var` changed with it),
and only the manifests' pinned gene-axis input differs between the #208 and #210
re-publishes. Each h5ad
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
| `disease-fusion/` | `bases.disease_deg.DiseaseFusion` (`catalog`, `manifest`, `diseases`, `record`, `deg2000`; digest-checked per read) |
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
`execution/make_gene_axis_docx.py` from the coverage CSV, the PNG, the six
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
python stage2/01_bases/execution/build_disease_fusion.py verify
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
CP10K for `G_2000`/`G_5000` from an interim pre-#190 axis. On the post-#194 axes the
measured tenv shares were **15.0548953% / 37.4273748%** (treated-only 15.0892505% /
37.3846128%); the ADR 0071 v2 axes re-measure **13.1626914% / 33.8379148%**
(treated-only 13.2072122% / 33.7895449%), which is what the current manifests record.
The ~2.5x denominator argument is unchanged in every generation and the manifests
record the released-axis values as authoritative. ADR 0066's figures are superseded,
not re-based.

**Rejected-source rationale (resolved by #212).** The `rejected[]` strings are
declaration data in `execution/disease_deg_sources.py`, copied verbatim into each
category's `catalog.json` and digest-pinned through the fusion, gene-axis and view
manifests. Two of the strings carried over from the pre-v2 registry were false in
the published products: `nHF`'s `GSE183852` entry claimed the CM labels live only
inside 2.1–12 GB Seurat `Robj` objects and that re-clustering is forbidden (the
authors' 6.7 MB `MOESM4` `Cardiomyocytes` sheet carries the labels, and the
measured v2 yield is 70 DEGs / 50 on-axis), and `nHF`'s third entry rejected
whole-tissue bulk for not being cardiomyocyte-derived (ADR 0071 tier 3 admits it).
`ICM`'s `GSE298023` entry shared the false re-clustering wording. #211 could only
correct the comment of record while the release was frozen; #212 lifted the
freeze, rewrote the three strings in the declaration, and re-published the cascade
(`disease-deg` → `disease-fusion` → `gene-axis` → both views). The catalogs now
carry the measured reasons; no other declaration byte, axis gene set, view payload
or rule-document byte changed (gate: [`log/acceptance.md`](log/acceptance.md) §12).
The item-local
`local/deg/*_deg.json` platform strings that #211 noted were **not** rebuilt here
(they are unrelated item artifacts, and rebuilding them would move every catalog's
pinned `derived_inputs` digests for no scientific gain); their author-log erratum
notes still say "rewritten on the next DEG rebuild".

**Cosmetic, unfixed (out of scope).** `bases/disease_deg.py` emits a pre-existing
`SyntaxWarning: invalid escape sequence '\'` from a #187 docstring; it is cosmetic and
owned by a follow-up.

## Accepted results and limits

- Replay is byte-identical for every released product (see `log/acceptance.md`); the
  two axes are **nested** (`G_2000 ⊂ G_5000`, 2,000 shared genes, Jaccard 0.4) and
  every axis index is in the frozen 18,164 space.
- `corpus/` carries no normalization, no delta and no HVG subsetting; `view-*` normalizes
  once on the **full-axis** library total and then subsets, so `X`/`control`/`Δ` share one
  scale across axes (ADR 0066).
- The axis is a selection on measurements, not a biological claim; `f2000 ∈ 0..6` is
  coarse, `rank_sum`'s scale is dominated by `AF`'s one-gene denominator (see the fill
  rule above), and the six diseases are not statistically independent (shared
  `GSE270788` and Reichart controls, declared not corrected).
- `02_models_tenv` still reads the retired `training-corpora/tenv.h5ad`; re-pointing it at
  `view-2k/`/`view-5k/` and re-running its accepted matrix are **deferred** (ADR 0066
  *Deferred*), as is the evaluation redesign and recomputing the tenv-HVG ranking from the
  corpus (the default `--training` already names that live payload).
- `GSE183852` remains the one closed `nHF` candidate: its authors' 6.7 MB `MOESM4`
  `Cardiomyocytes` sheet does carry the CM labels (so reading its `meta.data` would not be
  re-clustering), and the v2 definition measures 70 DEGs there, 50 on the frozen axis —
  2.5% of a 2,000-gene target, so the 8.1-12 GB `Robj` label extraction is not worth
  running. `nHF` therefore publishes its whole 1,661-gene union. Its padj-only yield
  (4,927 on-axis) would only matter under a softer fold tier, which v2 does not admit.

## Historical boundary

The pre-0066 measured truth, successor envelope, diagnostic PR chain and the
unweighted quota-free/direction-composition axis evidence are retired science: do not
resurrect their semantics. Their ADRs (0061/0062/0063/0065) and the research notes under
`docs/research/` keep their original bytes as history; this README describes only the
post-rebuild tree.

## Unified TOX independent admission candidate (ADR 0072)

Before TOX DEG analysis, run `execution/build_tox_admission.py publish` and
`verify` under the arc environment. This metadata-only build writes
`execution/checkpoints/tox-admission/{admission,manifest}.json`, preserving the
accepted disease release. It consumes `execution/tox_drug_evidence.json` and
`tox_source_admission.json`, pins their primary evidence bytes and chemical
annotation, and reads no new DEG result. `bases.tox_admission.ToxAdmission`
checks the candidate payload digest and optionally every input digest;
`drugs(source, admitted_only=True)` is the pooling-membership contract.
Names with conflicting structural identities fail closed. Independent author
negative labels, unresolved structures, antibodies/peptides and pending evidence
never enter TOX. A source marked admitted is eligible for tiered analysis, not
proof of a positive result. The ledger records finite discovery scope, excluded
platform duplicates and the LINCS full-transcriptome factual correction without
expanding the accepted source tiers.

The finite ADR 0072 source review admits twelve tier-2 designs, including the
reassessed antiviral source GSE247541 and the independently supported drug subset
of GSE262419. All panel identities remain inspectable; shared evidence transfers
only by full parent identity. DToxS official-label and primary-case snapshots are
declared under `pert_data/DToxS-iPSC/local/tox_evidence/`; the fetch recipe verifies
their pinned bytes. Case attribution uncertainty is retained. Pending panel
records carry the reviewed evidence-set pointer and date; missing support in this
finite set is neither a negative label nor a claim of exhaustive drug research.

### Unified TOX measurement candidate (#216)

`execution/build_tox.py publish|verify` computes only
`execution/checkpoints/unified-tox/tox/`; it does not replace the released
mechanism-class products. Run through `execution/capped.sh /usr/bin/time -v`
with the R configuration above. `verify` re-fits every entered dataset from the
pinned measurements into separate scratch space and compares all output and
manifest digests. Complete all sources in a tier before testing the cumulative
axis-resolved 2,000-DEG stop rule.

`bases.tox.ToxCandidate` inspects the candidate `DEG2000` and source membership
audit; the ordinary `DiseaseDeg` reader reads the candidate root as well. There
is one DEG/rank pair per independent dataset, a full independent admission
snapshot, exact original treated/vehicle/excluded well audit, and one TOX
`tox.deg2000.csv`. Mechanism is annotation and never introduces extra source
votes. Original DToxS edgeR, GSE243674 donor-paired DESeq2 and ScreenSeq Welch
tests are retained. Additional count sources use the shared limma-voom worker;
FPKM is explicitly log2(FPKM+1) on the existing expression-limma path, while
fractional estimated counts use the existing documented integerisation helper
with deviation/scale statistics. Runtime peaks remain in the work directory,
outside deterministic scientific payloads. Integrated publication, axes and
views follow only after successor acceptance.
