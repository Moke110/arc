# Stage 2 shared foundations

> **Superseded layout, rebuild pending (ADR 0066).** This step's `output/` tree is being
> rebuilt as `raw_counts/` → `corpus/` → `view-2k/`/`view-5k/`, with `disease-deg/` and the
> released `gene_axis_2000.csv` / `gene_axis_5000.csv` in `gene-axis/`. The measured truth,
> the successor envelope, the diagnostic PR figure and `log/` are retired. Everything below
> describes the pre-0066 tree until the rebuild lands.

## Rebuild foundation: `raw_counts/` and `corpus/` (ADR 0066, #186)

The rebuild's first slice has landed. These two product areas are the new
`source/` -> axis chain; the pre-0066 products above stay in place read-only until
the acceptance/removal ticket (#192) deletes them.

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
source-accepted tenv/Tahoe payloads, re-checks their bytes against
`log/records/native_tenv_verify.json` and
`log/records/tahoe_counts_source_acceptance.json`, re-derives the Tahoe token
mapping from the immutable source `gene_metadata.parquet`, and records
`provenance.mode: source-accepted-replay` with a pointer to the acceptance record
and the recorded source-replay evidence. The stronger mode
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
re-exports these under the old names so the accepted pre-0066 evidence keeps
working.

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
n_measured, test_method, member_drugs, n_small_molecule_drugs` for the fusion
ticket.

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
never treated as independent experiments.

```bash
python stage2/01_bases/execution/build_disease_deg_cardiotox.py publish
python stage2/01_bases/execution/build_disease_deg_cardiotox.py verify
```

The build uses the repository's R build environment
`/home/chang/miniconda3/envs/biomni_e1/bin/Rscript` (R 4.5.3, edgeR 4.8.2,
DESeq2 1.50.2; recorded as `conda:biomni_e1`), because the `arc` env's R has
neither package. The pinned DToxS pipeline calls `decideTestsDGE`, removed in
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

## Layout and the single public implementation

| path | role |
|---|---|
| `output/api/bases/` | the single public reader/metric implementation (tracked code) |
| `output/gene-identity/` | indexed identity, `genes.csv` and the frozen `tables/` |
| `output/tenv-raw-counts/` | tenv integer raw counts, probe evidence and seal |
| `output/tahoe-raw-counts/` | projected Tahoe integer raw counts, token evidence and seal |
| `output/shared-expression/` | verified CP10K/log1p shared projection (reconstruction intermediate) |
| `output/training-corpora/` | current separate Tahoe and tenv training releases |
| `output/historical-training-release/` | historical v1 truth-fit release — tenv h5ad, seal and exclusions; the rejected Tahoe payload was removed by #157 |
| `output/tenv-truth/` | accepted measured 1377-contrast tenv truth |
| `output/tenv-truth-successor/` | accepted exact successor envelope for the current training release |
| `output/tenv-hvg/` | descriptive and three train-only HVG rankings plus partitions |
| `figure/diagnostic_tenv_pr.png` + `figure/data/…csv` + `figure/script/…py` | report panel, complete plotted values and the CSV-only renderer |
| `execution/` | build, audit, migration, figure-provenance and test code; disposable state in `execution/checkpoints/` |
| `log/records/`, `log/health/`, `log/diagnostic/`, `log/historical/` | acceptance/execution, repo-health, dissolved and historical evidence |

Add `output/api` to the import path and use the published modules; execution
commands do this through `execution/public_api.py`:

```python
import sys; sys.path.insert(0, 'stage2/01_bases/output/api')
from bases.identity import GeneIdentity
```

## Identity

`output/gene-identity/genes.csv` fixes `gene_index` 0–18163 in the original
four-CM-source order. `output/gene-identity/tables/genes.csv` retains the exact
original axis bytes; the crosswalk retains source evidence. FAM153B is index
5091: tenv ENSG00000182230 and Tahoe ENSG00000289731/token61482.

`indices(ids, source=...)` resolves native namespaces with explicit ambiguity
errors; `labels` and `native_ids` cross the human/native boundary. The shared
mapping files are derived solely from raw gene metadata and pinned HGNC evidence.
The sealed identity product still binds the historical v1 mapping digests, which
are preserved at `log/historical/shared-v1-mappings/` (see Limitations).

```bash
python stage2/01_bases/execution/build_identity.py
python stage2/01_bases/execution/build_identity.py --verify
```

## Corpora: axes and sample counts

Two separate training corpora are retained on the same frozen 18,164-gene axis,
and both hold **treated rows only** — vehicle controls are not rows; they enter
through `layers/control`. `tenv` is a single-cell-line hiPSC-CM TempO-Seq corpus;
`tahoe_bulk` is a 50-cell-line pseudobulk corpus derived from the Tahoe-100M
atlas. Both are single-time-point experiments (tenv records `time_h = 24`; the
retained Tahoe metadata carries no time axis), so time is not one of the counted
axes.

Axis coverage of `output/training-corpora/`:

| corpus | 细胞 (cells) N_c | 药物/细胞 (drugs per cell) N_drug | 剂量/(细胞，药物) (doses per cell×drug) N_dose | 重复数/(细胞，药物，剂量) (replicates per cell×drug×dose) N_rep | 总样本量 (total samples) N = N_c·N_drug·N_dose·N_rep |
|---|---:|---:|---:|---:|---:|
| `tenv` | 1 (hiPSC-CM) | 459 | 3 (0.2 / 1 / 10 µM) | 3.05 mean (mode 3; 1,320×3, 48×4, 9×6) | **4,206** |
| `tahoe_bulk` | 50 (cell lines) | 352 | 3 (0.05 / 0.5 / 5 µM) | 1.24 mean (mode 1; 43,036×1, 7,941×2, 1,400×3, 350×4, 50×16) | **65,318** |

`N` is the realized count: the sum, over realized (cell, parent, dose) strata,
of the number of training profiles in each stratum. The schema
`N = N_c · N_drug · N_dose · N_rep` is therefore a full-grid estimate, not an
identity — `N_c · N_drug · N_dose` counts full-grid combinations and `N_rep` (mean) is the
average number of profiles per realized stratum. Both corpora are imbalanced:

- **`tenv`** — all `1 × 459 × 3 = 1,377` strata are realized; the three replicate
  wells give a 4,131-profile grid, raised to 4,206 by 48 fourth wells and 9
  six-profile strata (two source compound names share one parent: Tamoxifen /
  Tamoxifen citrate, ammonium perfluorooctanoate / perfluorooctanoic acid,
  potassium perfluorooctanesulfonate / perfluorooctanesulfonic acid).
- **`tahoe_bulk`** — the full grid is `50 × 352 × 3 = 52,800` combinations; 52,777 are
  realized because 23 (cell, parent) pairs were profiled at two doses, not three.
  The 65,318 profiles are `43,036×1 + 7,941×2 + 1,400×3 + 350×4 + 50×16`; the
  multi-profile strata come from the 24 parents that aggregate 2–3 source drug
  names and from conditions replicated across plates (Adagrasib at 0.05 µM has
  16 wells, the only size-16 stratum).

How each axis is defined:

- **Cell** — `obs.cell_identity`. tenv has one cell identity (`hiPSC-CM`;
  source `cell_type = ipsc_cm`); Tahoe is the 50 atlas cell lines
  (`obs.cell_name`, e.g. `A549`).
- **Drug** — the desalted parent structure (RDKit `Cleanup → FragmentParent →
  Uncharger`, stereochemistry retained), i.e. one `drug_parent` value. Source
  names collapse into parents: tenv keeps 462 source compound names → 460 source
  `drug_key`s → **459 parents**; Tahoe keeps 377 source drug names (379 treated
  names before structure exclusion) → **352 parents**.
- **Dose** — the delivered concentration in µM (`obs.dose_um`), exactly the three
  values listed per corpus.
- **Replicates** — the number of training profiles per realized (cell, parent,
  dose) stratum. This counts replicate wells and salt/alias forms of the same
  parent; it is not a pure technical-replicate count.

Row populations behind the table (controls and blanks are counted separately and
never as drugs):

| corpus | native raw profiles | shared-expression rows | training rows | vehicle controls (not rows) | excluded |
|---|---:|---:|---:|---:|---|
| `tenv` | 4,510 (16 plates) | 4,478 | **4,206** | 272 (192 same-plate DMSO matched; 80 Media excluded) | 32 technical blanks |
| `tahoe_bulk` | 67,018 (14 plates) | 67,018 | **65,318** | 1,400 `DMSO_TF` (all matched) | 300 missing/invalid structure |

## Expression preprocessing: from `pert_data/` raw counts

Every value the models read is a deterministic function of integer counts under
`pert_data/`; no source file is modified, and every value in the training h5ads
can be recomputed from the raw counts. The chain is raw integer counts → native
raw-count matrix → frozen-axis CP10K/log1p projection → training row selection,
matched-control mean and delta.

| stage | tenv | tahoe |
|---|---|---|
| 0 raw counts (`pert_data/…/source`) | 16 TempO-Seq plate CSVs + well hash table; 4,558 wells | 3,388 sparse parquet shards; 95,624,334 cells; 14 plates × 96 wells |
| 1 native raw counts (`output/*-raw-counts/`) | 4,510 × 19,613 `int32` | 67,018 × 18,164 `int32` |
| 2 shared expression (`output/shared-expression/`) | 4,478 × 18,164 `float32` | 67,018 × 18,164 `float32` |
| 3 training corpora (`output/training-corpora/`) | 4,206 × 18,164 `float32` × 3 layers | 65,318 × 18,164 `float32` × 3 layers |

### 0. Raw counts in `pert_data/` (immutable)

- **tenv** — `pert_data/GSE262419/source/` (GEO GSE262419, Tsai et al. TempO-Seq
  hiPSC-CM environmental-chemical panel): 16 gzipped per-plate matrices
  `GSE262419_Plate01.csv.gz` … `Plate16.csv.gz`, one column per well and a
  `Genes` column of TempO-Seq probe identities; `GSE262419_hash.csv.gz` maps each
  (plate, well) to a chemical name, concentration in µM and solvent (`VEH` =
  DMSO, `Media`, `Blank`). Values are integer counts per probe and well. The
  hash table declares 4,558 wells; the 48 `ET-1` wells (peptide, 0.01 µM) are
  outside the admitted small-molecule set, leaving the 4,510 retained wells.
- **tahoe** — `pert_data/Tahoe-100M/source/` (HuggingFace `tahoebio/Tahoe-100M`,
  rev `349e2a55…`, CC0-1.0): 3,388 `data/train-*.parquet` shards holding one row
  per cell with sparse `genes` (token ids; the first token of each list is a
  marker) and `expressions` (integer UMI counts); metadata tables
  `sample_metadata.parquet` (1,344 wells), `cell_line_metadata.parquet` (50
  lines), `gene_metadata.parquet` (token ↔ Ensembl) and `drug_metadata.parquet`;
  `manifest.sha256` seals the bytes. The 14 plates × 96 wells contain 28
  `DMSO_TF` vehicle wells and 1,316 treated wells.

Neither source stores CP10K, log values, pseudobulk or matched-control means;
`pert_data/` is raw counts only. Both items carry the external-data provenance
contract (`author_log.json` plus a `manifest.sha256` over the downloaded bytes).

### 1. Native raw-count matrices

- **tenv** (`execution/build_native.py`, `bases/native.py`) — every probe is
  resolved to a pinned HGNC gene using `pert_data/Annotation/source/hgnc_complete_set.txt`
  (the declared local `pert_data/Annotation/local/cm_compound_annotation.csv` supplies
  the per-compound structure and pharmacology annotation): 22,537 probes, of which 21,733
  match directly, 705 by previous symbol, 17 by Entrez id, 3 by alias and 79 are
  unmapped. Probes that resolve to the same HGNC are **summed**, giving 19,613
  native genes; `X` is `int32` counts for 4,510 wells (4,206 treated / 272
  control / 32 blank) with `lib_size` and `n_genes_detected`; probe evidence is
  `output/tenv-raw-counts/tenv_probe_mapping.csv`. No normalization.
- **tahoe** (`execution/prepare_tahoe.py`, `bases/tahoe.py`) — each cell's
  non-marker tokens are mapped through `gene_metadata.parquet` and the sealed
  crosswalk to the frozen axis, and integer UMI counts are **summed** per
  (plate, well, cell line) into one pseudobulk profile; marker and off-axis
  tokens are dropped. The result is 67,018 profiles × 18,164 genes `int32`
  (65,618 treated + 1,400 `DMSO_TF` controls), with `n_cells` per profile and
  token evidence in `output/tahoe-raw-counts/token_mapping.csv`. Tahoe is
  projected at this stage; its native 62,710-token vocabulary is not retained.

Both raw releases are verified by bounded replays: every tenv probe and every
one of the 95,624,334 Tahoe cells.

### 2. Frozen-axis projection, CP10K and log1p (shared expression)

`execution/build_shared.py` / `bases/raw.py` project every raw profile onto the
ordered 18,164-gene axis and normalize each profile independently:

```
X[p, g] = ln(1 + 10^4 · C[p, g] / Σ_{g' ∈ axis} C[p, g'])
```

- The axis is `output/gene-identity/tables/genes.csv` (sha256 `19ce44ac…`):
  18,164 rows of `symbol, hgnc_id, entrez_id, ensembl_id`; `match_type` is direct
  17,067 / alias 1,086 / entrez 11.
- Mapping: tenv resolves all 18,164 genes by exact HGNC and drops its 1,449
  native off-axis genes; Tahoe resolves 18,162 by exact Ensembl, 1 by
  corroborated Ensembl and 1 through the crosswalk (FAM153B, `gene_index` 5091:
  tenv `ENSG00000182230`, Tahoe `ENSG00000289731` = token 61,482).
- The library total is the sum over exactly the 18,164 projected genes, never
  the native full-transcriptome total. `ln` is the natural log (`np.log1p`);
  values are stored as `float32`; a zero count stays exactly `0.0`.
- Rows: tenv 4,478 (4,206 treated + 272 controls; 32 blanks dropped), Tahoe
  67,018 (65,618 + 1,400). No library-size or detected-gene filter is applied.

Independent replay recomputes all 4,478 and 67,018 rows from raw counts with
`max_abs_error` 0.0 ([`log/records/shared_v2_verify.json`](log/records/shared_v2_verify.json)).

### 3. Training corpora: rows, matched controls, `X` / `control` / `delta`

`execution/build_corpora.py` / `bases/corpora.py` keep treated profiles only and
attach a matched-control mean:

- **Rows** — tenv 4,206; Tahoe 65,318 (65,618 treated profiles minus 300 whose
  structure is missing or invalid: `Sacubitril/Valsartan` and `Verteporfin`, each
  at 3 doses × 50 cell lines). Exclusions are recorded in
  `output/training-corpora/*_exclusions.csv`.
- **Matched controls** — tenv: same plate, `DMSO` only (192 wells, 12 per
  plate); Media controls are not eligible. Tahoe: same plate and same cell line,
  `DMSO_TF` (1,400 profiles).
- **`X`** — the treated row of the shared expression matrix.
- **`layers/control`** — arithmetic mean of the matched controls' rows:
  `control[s] = mean_{c ∈ C_s} X[c]`. Each control is first normalized on its own
  projected library, then averaged in `float64` and stored as `float32`; control
  libraries are never pooled before normalization.
- **`layers/delta`** — `X − control`, stored as `float32`.

The h5ad `var` is the frozen axis plus `gene_index`; `obs` carries the source
profile/well identity, `plate`, cell identity, source drug name and parent
InChIKey, `dose_um`, `source_inchikey`, `source_parent_key`, matched-control ids
and a `group_id` (parent–dose hash). Training files contain no significance
labels.

Deliberately absent from every layer: batch or plate correction beyond the
matched-control subtraction, per-gene centering or scaling, HVG/gene subsetting
(all 18,164 genes are kept), quantile or median normalization, imputation,
control rows, and any delta in the raw or shared products.

Every `X`, `control` and `delta` value is replayed from raw counts with
`max_abs_error` 0.0: 4,206 rows / 192 controls for tenv and 65,318 rows / 1,400
controls for Tahoe
([`log/records/training_v2_verify.json`](log/records/training_v2_verify.json)).

## Native raw counts

`build_native.py` streams the 16 GSE262419 plate files one plate at a time, sums
duplicate probes by HGNC, and writes `output/tenv-raw-counts/tenv.h5ad` with 4510
raw wells and 19613 native genes (digest `efd07632…`). It retains controls and
blanks and has no normalized expression or significance layers. The manifest
seals every immutable plate, hash table, HGNC table, compound registry and output
file. Independent verification sums the source probes directly for every raw well.

`output/tahoe-raw-counts/tahoe.h5ad` holds projected 18164-gene Tahoe integer
pseudobulk counts for 67018 profiles; `token_mapping.csv` keeps the original
token evidence and the manifest seals all 3388 source shards. Independent source
acceptance replayed every shard, profile and cell (95624334).

```bash
python stage2/01_bases/execution/build_native.py --verify
python stage2/01_bases/execution/prepare_tahoe.py verify
```

The CLIs accept `--data-root`, `--output`, `--tenv`, `--tahoe` and `--tables` for
independent worktrees. `bases.raw.tahoe_axis_mapping(axis, crosswalk,
gene_metadata)` exposes the original token-to-index evidence without loading
expression. Tahoe's raw preparation boundary projects directly onto the frozen
18164 genes; it does not require recreating the deleted 62710-gene matrix.

## Shared expression and training corpora

`output/shared-expression/` holds the bounded CP10K/log1p projection
(67018 × 18164 and 4478 × 18164 float32); independent verification replays every
projected row with `max_abs_error` 0.0. `output/training-corpora/{tenv,tahoe_bulk}.h5ad`
are the current separate training releases and keep individual treated wells.
`X` is natural log1p of raw counts projected onto the full fixed axis and scaled
to 10000. `layers/control` averages independently normalized controls;
`layers/delta` is treated minus control. Tenv uses same-plate DMSO, excludes
Media controls and has cell identity `hiPSC-CM`. Tahoe uses same-plate,
same-cell-line DMSO_TF. Training files contain no truth labels.

Samples retain stable source profile IDs, source salt structures/keys,
stereochemistry-aware normalized parent keys, micromolar doses and exact
matched-control IDs. Parent identity uses RDKit Cleanup, FragmentParent and
Uncharger. The tenv release contains 4206 wells and 1377 parent–dose groups;
81 source-key conversions are explicit. Blank and missing-structure exclusions
are recorded. Group IDs are independent of partition membership.

```bash
python stage2/01_bases/execution/verify_training_source.py --source tenv
python stage2/01_bases/execution/verify_training_source.py --source tahoe
```

Replay joins immutable profile/control IDs, so rebuilding a raw matrix in a
different row order cannot change identity. Every treated/control/delta value is
recomputed from raw counts in bounded windows. Equivalent tautomer SMILES
serializations are audited with mandatory exact InChIKey and stereo-aware
parent-key equality; original training strings remain unchanged. New raw files
have new byte seals and are not claimed byte-identical to deleted files.

`build_corpora.py` builds the pair after both bases-native inputs and the
bounded shared expression release exist. `bases.raw.build`/`verify` own the
shared projection/normalization; explicit source subsets support incremental
preparation. All defaults point inside this step. `training_header` opens small
h5ad metadata only, never truth files.

## Measured truth and the exact successor binding

`build_significance.py --native <tenv raw counts>` builds separate parent–dose
truth from raw counts on exactly 18164 genes; the product is `output/tenv-truth/`
(1377 contrasts, 4206 treated wells, 273 tested-empty, 6 multi-plate). The
bases-owned DESeq2 worker uses a two-sided Wald test, plate effects for
multi-plate contrasts and within-contrast BH. The partition-authorized reader
exposes only requested groups. Untestable and tested-empty are different states.
Truth fitting and full dual-corpus source acceptance have their own records; a
successful training build alone does not establish either result.

`output/tenv-truth-successor/` is the accepted exact-equivalence envelope that
authorizes the current training release with those labels. It is re-issued after
a path move, never edited: the superseded envelope bytes are preserved at
`log/historical/seals/truth-successor-envelope-v1/` and the new digest is
recorded in [`log/156-successor-reissue.json`](log/156-successor-reissue.json).
The proof compares integer contrast counts, sample/control identity, the
statistical recipe and bit-exact `X`/`control`/`delta` layers; it performs no new
fit.

Run heavy commands serially under the repository memory caps
(`execution/capped.sh`). Tests are individual `unittest discover` files under
`execution/tests/`.

## HVG, metrics and the report figure

HVG selection, retrieval metrics and evaluator commands are documented in
[`execution/hvg_metrics.md`](execution/hvg_metrics.md). Model selection uses
training-only rankings from `output/tenv-hvg/`; full-tenv descriptive ranks are
not model inputs.

The report panel is `figure/diagnostic_tenv_pr.png`. Its complete plotted values
are `figure/data/diagnostic_tenv_pr.csv` — curve coordinates, fixed
contributing/eligible group counts, the displayed equal-group macro AP per panel
and fold, the fixed-group definition, and the uncertainty/reference statement.
`figure/script/render_diagnostic_tenv_pr.py` reads only that CSV and reproduces
the PNG byte for byte (sha256 `35d14266…`); it never opens `output/`.

The panel covers the three test folds, 4206 measured wells and 1377 parent–dose
groups and is measurement-derived diagnostic QA, not trained-model performance.

## Gene-axis decision evidence (#183)

`output/gene-axis/` is the **decision-preparation** evidence for the frozen
18,164-gene axis: per-source HVG and DEG ranks, an unweighted rank-sum fusion,
source overlap, direction composition, an explicit inflection answer and the
tenv non-HVG curve. It chooses nothing — the axis size `N`, the per-source quota,
the direction rule and the animal-source weight remain #173's decision. It builds
no corpus, trains no model, publishes no `genes.csv`, and never opens a truth or
Heliste label file; the only tenv input is the published descriptive HVG ranking.

Public implementation: `output/api/bases/gene_axis.py` (pure functions: symbol /
alias / Ensembl alignment, Seurat-style normalized dispersion, the continuous DEG
effect score, within-source ranks, rank-sum fusion, overlap, direction
composition, Kneedle + second-difference inflection, tenv non-HVG curve, and the
gene-index helpers that keep the sorted fused table from being read by row
position). Execution: `execution/build_gene_axis.py` (build + `--verify` replay),
`execution/gene_axis_sources.py` (the two `pert_data/` cardiotoxicity adapters)
and `execution/build_mouse_orthologs.R` (the pinned mouse→human table; it writes a
candidate and replaces the accepted table only after a sha256 check, with
`--args --accept` for a deliberate re-pin).

Fused sources (one rank column each, tenv included): the published tenv HVG rank,
GSE141910 and GSE116250 (HF), GSE243674 (hiPSC-CM anthracyclines), GSE128188 (AF),
GSE282618 (hypertensive heart disease, exploratory), the two mouse MIRI DEG-only
items GSE160516/GSE108940 (homologene-mapped), and the two `pert_data/`
cardiotoxicity sources DToxS-iPSC (22 upstream tables) and GSE244740 (ScreenSeq
counts re-derived). `GSE137293` is listed as supplementary-only and `LINCS-LDS-1593`
as a coverage note; neither is fused.

Declared recipe (details and constants in `output/gene-axis/manifest.json`):
per-source HVG = normalized dispersion on the item's own deposited log expression
over its **whole measured transcriptome** (z-scored in 20 equal-width
`log1p(mean)` bins, singleton bins score 1, undefined sorts last, source-row
tie-break — the `bases/hvg.py` family without `expm1` recovery because the disease
matrices are heterogeneous log scales); a source's declared HVG set is its top-K
source rows, and `hvg_share = n_hvg_on_axis / n_hvg_total` (share of the source's
HVG set that lands on the axis), with `n_hvg`, `n_hvg_total`, `n_hvg_on_axis` and
`hvg_density = n_hvg_on_axis / n_axis_expression` reported; DEG =
`|log2FC| * -log10(padj)` with NA padj treated as non-significant and multi-contrast
max-collapse; unmeasured genes get rank `n_measured + 1` and rank-normalized value
1.0; fusion = unweighted sum of `rank / n_measured` (the animal half-weight column
is a reported sensitivity, not a chosen rule); fused order ties break on
`gene_index`. Top-N direction composition is reported under two **different**
orderings: `symmetric` (ascending fused rank) and `max_abs_log2fc` (descending
max cross-source |log2FC|, fused-rank tie-break). Multi-way overlap intersections
report the size and true fused ranks with an explanatory `note`; Jaccard is left
blank because the operands are identical.

Headline evidence (exact digests in the manifest, numbers in
[`log/gene-axis/record.md`](log/gene-axis/record.md)):

- Axis coverage 0.68–1.00 per source; per-source HVG share (of the source's own
  top-2,000 HVG set, on the axis) 0.1465–0.9375 for the expression sources (tenv
  1.0 by definition, density 0.1101) with `hvg_density` 0.0162–0.1389; DEG
  up/down shares range from 0.000/0.000 (GSE108940: one significant gene at
  padj<0.05) to 0.53/0.45 (DToxS).
- **No inflection.** Kneedle (online `S=1.0`) puts its broad arc maximum at
  `N≈1314`; the second-difference peak sits in the tail (`N≈18156`); the criteria
  do not agree, so `inflection.json` reports `"inflection": false` and `n: null`.
  All ten leave-one-source-out replays also report no inflection.
- Tenv non-HVG share of the **true fused top-N** at declared `K=1000 / 2000 /
  5000`: 0.884 / 0.782 / 0.540 at N=500 and 0.913 / 0.823 / 0.610 at N=2,000,
  rising to the axis background 0.945 / 0.890 / 0.725 at N=18,164. The fused top
  is enriched for tenv HVGs relative to the background but remains mostly
  (82.3% at K=2000) outside tenv's top-2,000.
- HF sources agree (GSE141910 ∩ GSE116250 Jaccard 0.4255, 5,320 genes); the
  AF/HHD chronic-remodelling pair is nearly disjoint (0.0043) as are HF×remodelling
  pairs (cross-pair mean 0.0233).
- Top-N direction composition is 52–54% up / 46–48% down under symmetric
  inclusion; the max-|log2FC| ranking is 71–88% up at the top. Across all
  significant genes the down set has the better median fused rank (8,031 vs
  8,779 for up).
- Effect×replication rank vs DEG-count rank: Spearman 0.81 (ties broken by
  ascending fused rank) and top-N overlap 0.62 / 0.56 / 0.48 / 0.45 / 0.76 at
  N = 500 / 1000 / 2000 / 3000 / 5000.

```bash
python stage2/01_bases/execution/build_gene_axis.py
python stage2/01_bases/execution/build_gene_axis.py --verify
R_HOME=/home/chang/miniconda3/envs/arc/lib/R   /home/chang/miniconda3/envs/arc/lib/R/bin/exec/R --no-save --slave   --file=stage2/01_bases/execution/build_mouse_orthologs.R   # only to regenerate the pinned table
python stage2/01_bases/figure/script/render_gene_axis_source_summary.py   --csv stage2/01_bases/figure/data/gene_axis_source_summary.csv   --png stage2/01_bases/figure/gene_axis_source_summary.png
```

The figure is `figure/gene_axis_source_summary.png` (four panels: per-source
HVG/DEG shares, fused rank-sum histogram with the candidate inflection, tenv
non-HVG share vs N, source-overlap Jaccard matrix); its complete plotted values
are `figure/data/gene_axis_source_summary.csv` and the CSV-only renderer
reproduces the PNG byte for byte (sha256 `91673888…`) without opening `output/`.
`execution/tests/test_gene_axis.py` is fixture-driven and includes the
truth-path poisoning test plus index-space regression fixtures where fused-rank
order differs from `gene_index` order.

## Shared `@d` evaluation and figure APIs

`output/api/bases/evaluation.py` owns the ADR 0064 reliability family used by the
tenv comparison; `02_models_tenv` imports it and re-implements none of it. The
direct family is unchanged.

- **Label side.** `measured_effect_gate(padj, log2_fold_change)` marks every gene
  with `padj < 0.05`, direction `sign(log2FC)`, and returns
  `MeasuredEffects(signed, magnitude, indices)`. The declared constants are
  `PADJ_THRESHOLD`, `THRESHOLD_GRID_POINTS = 100`, `MIN_HITS_FOR_SPEARMAN = 3`,
  `SPEARMAN_NEUTRAL = 1.0`, `GENE_COUNT` and `RELIABILITY_METRICS`.
- **Predicted side.** `predicted_rank_table(delta_mean)` is the full-axis
  `|Δ̂_mean|` rank table, retained unchanged from the previous comparison.
- **Scan.** `threshold_grid(M)` returns `d = k · 0.01 · M` for `k = 1..100`;
  window membership is `|Δ̂_mean| >= d`, so the strictest point keeps the window
  of the gene that defines `M` non-empty instead of emptying every window. It
  does not select every group's strongest gene: at `d = M` only the group(s)
  whose maximum is `M` can be non-empty (12/459 fold 0, 33/459 fold 1 here), so
  the top of the grid is a degenerate edge.
- **Metric functions.** `precision_at_d`, `recall_at_d`, `f1_at_d` (harmonic),
  `spearman_at_d` (neutral `mult = 1` below three hits) and `reliability_at_d`
  are pure functions over counts and arrays. `reliability_curve` binds one
  group's window to the grid; `aggregate_reliability` is the equal-weight macro
  over a fold's test groups with the null/changed split. The family is **total**:
  an empty window on a null group scores `1`, an empty window on a group with
  measured effects scores `0`, a non-empty window on a null group scores `0`.
  No metric returns NaN and none carries a defined-group count.

`output/api/bases/plots.py` owns the complete plot-value schema and the drawing
API both step figures use. New panel kinds added for ADR 0064 are a `curve`
panel with a per-run threshold marker (`d*`) and an optional reference line
(`F1@d` on the Reliability panel), and a three-class `heatmap` panel; the `bar`
panel is reused for the direct metrics and the reliability summaries. The
figure row records its panel grid, so a replay is driven by the recorded layout
rather than a hard-coded shape. `output/api/bases/palette.py` gains the
`reference` and `marker` roles and a recorded sequential heat-map ramp; entry
hues, line/fill shades and the learning-curve training/validation families are
unchanged. `02_models_tenv/figure/script/render_figures.py` is the only step
consumer.

The retired k-parameterized indirect family and its undefined-group machinery
are no longer called by any live step; `bases/metrics.py` (the earlier retrieval
module) remains untouched and is not imported here.

## Accepted results and limitations

Results reused from the accepted records (#156 re-verified digests and path
bindings rather than repeating the scientific computation, except where a bounded
gate is named in [`log/156-gates.md`](log/156-gates.md); #157 completed the final
audit, the removals and a fresh verification pass recorded in
[`log/157-final-summary.md`](log/157-final-summary.md)):

- Gene identity: 18164 genes, one frozen axis, FAM153B at index 5091; full-axis
  HGNC round trip and both source namespaces verified.
- Raw counts: tenv 4510 wells (4206 treated / 272 control / 32 blank); Tahoe
  67018 profiles over 14 plates with 1400 controls; 3388 source shards accepted.
- Shared expression: 67018 and 4478 rows replayed with zero numerical error.
- Training: 65318 Tahoe and 4206 tenv samples in the current release; both
  source replays pass with `max_abs_error` 0.0.
- Truth: 1377 eligible contrasts on the exact 18164-gene BH family, fitted from
  the historical v1 tenv release; three independent raw R refits accepted.
- Successor binding: full 1377-contrast / 4206-sample exact equivalence proof,
  re-issued against the stable paths.
- HVG: four complete 18164-gene rankings; 4206 wells split into three
  group-isolated folds (train sizes 2508/2508/2535) with the held-out isolation
  intervention re-checked.
- Diagnostic panel equal-group macro AP (combined / upregulated / downregulated):
  fold 0 0.0798 / 0.3349 / 0.1374, fold 1 0.0878 / 0.3674 / 0.1406,
  fold 2 0.0929 / 0.3522 / 0.1452.

Limitations and boundaries:

- The panel's displayed AP is the exact equal-group macro AP, not the area under
  the sampled fixed-group pooled curve. No model predictions survive in this
  scope, so these numbers are QA of the measured delta layer.
- Tahoe has no significance labels or PR and publishes none.
- Historical v1 **Tahoe** is rejected and superseded, not an accepted training
  input: its source replay misattributed and dropped counts. #157 removed its
  11.20 GB payload after proving no accepted product, envelope or live consumer
  binds it; the historical manifest still declares its digest
  `6e3ac2ec…`, and the rejection records and compact counterexample remain in
  [`log/records/`](log/records/) ([`log/157-removals.json`](log/157-removals.json)).
- The historical v1 **tenv** release remains the original truth-fit input; the
  successor envelope, not a CLI default, authorizes current labels.
- The sealed identity product keeps the historical v1 mapping vintage: 229 Tahoe
  `native_symbol` strings differ from the current annotation, with no axis
  change. Re-sealing is a deferred user decision (#156 D2).
- `historical-training-release/manifest.json → inputs.shared_manifest`
  (`8ac09bcb…`) seals a shared-v1 manifest that no longer exists; it is a
  documented historical gap, not a product defect.
- No stochastic tolerance is claimed: the diagnostic and successor replays are
  deterministic, and file-integrity hashes are recorded separately from metric
  values.

## Rebuild from the retained immutable inputs

The raw inputs have separate owners: tenv is
`output/tenv-raw-counts/tenv.h5ad`, Tahoe is `output/tahoe-raw-counts/tahoe.h5ad`,
and the frozen axis/crosswalk live under `output/gene-identity/tables/`. `--tenv`,
`--tahoe` and `--tables` override each binding independently. The legacy
`--native DIR` option selects both h5ads plus `DIR/tables`; it is useful for
small fixtures. Files are read in place, without copying large raw payloads into
a common directory.

Build commands refuse existing destinations. Run each heavy command under the
repository memory/thread caps, one at a time, after obtaining a resource slot.
These commands form the supported dependency order, not evidence that a pending
full build has passed:

```bash
python stage2/01_bases/execution/build_native.py
python stage2/01_bases/execution/build_native.py --verify
python stage2/01_bases/execution/prepare_tahoe.py build
python stage2/01_bases/execution/prepare_tahoe.py verify
python stage2/01_bases/execution/build_shared.py
python stage2/01_bases/execution/build_shared.py --verify
python stage2/01_bases/execution/build_corpora.py
python stage2/01_bases/execution/build_corpora.py --verify
python stage2/01_bases/execution/verify_training_source.py --source tenv
python stage2/01_bases/execution/verify_training_source.py --source tahoe
```

Existing accepted releases can be verified and reused rather than rebuilt. For
an independent worktree, use explicit canonical paths for raw sources, tables,
shared release, training release and drug metadata: defaults resolve to that
checkout, which does not contain the large immutable inputs. The identity release
uses the surviving `shared-v1` mapping CSVs, independent of their absent
expression payloads; `--shared <step>/output/shared-expression` selects the
current mappings when a re-seal is deliberately requested.

## Historical boundary

The historical `output/historical-training-release/manifest.json` seals its
original build. It does not describe a current shared-expression payload, which
was deleted. Separate replay records bind raw SHA, axis SHA, training SHA and the
selected training manifest without rewriting historical bytes. The current Tahoe
training release passes raw-source replay for all 65318 treated profiles and 1400
matched controls with zero numerical error.

The two-cell counterexample and actual v1/v2/raw evidence are in
[`log/records/tahoe_historical_raw_diagnosis.json`](log/records/tahoe_historical_raw_diagnosis.json),
with the failed replay in
[`log/records/tahoe_training_v1_rejection.json`](log/records/tahoe_training_v1_rejection.json).
The v1 Tahoe payload is not repaired, relabeled or claimed equivalent to the
current release; #157 removed its bytes after the dependency proof, keeping the
declared digest, the sealed manifest and that counterexample
([`log/157-removals.json`](log/157-removals.json) R1). The v1 **tenv** h5ad, which
the truth fit does depend on, is untouched.

The exact successor envelope, partition reader and acceptance commands are
documented in [`execution/truth_successor.md`](execution/truth_successor.md).

## Records, migration and rollback

[`log/`](log/README.md) holds the reproducible audit records. The migration that
produced this layout is the committed, idempotent
[`execution/migrate_products.py`](execution/migrate_products.py); its record is
[`log/156-migration-record.json`](log/156-migration-record.json), including the
old→new pairs, per-file inode/digest evidence and rollback steps:

```bash
python3 stage2/01_bases/execution/migrate_products.py --rollback \
    stage2/01_bases/log/156-migration-record.json
```

The metadata binding audit is
[`execution/audit_metadata_bindings.py`](execution/audit_metadata_bindings.py)
(successor of the #155 record script); the one-off payload re-hash is
[`execution/verify_payload_digests.py`](execution/verify_payload_digests.py).

#157 closed the migration: [`log/157-final-summary.md`](log/157-final-summary.md)
is the final scientific summary, [`log/157-gates.md`](log/157-gates.md) the gate
and health record, [`log/157-removals.json`](log/157-removals.json) every removed
payload with its digest and resolution proof, and
[`log/157-evidence-ledger.json`](log/157-evidence-ledger.json) the reused-vs-new
ledger. It also regenerated
[`log/156-metadata-binding-audit.json`](log/156-metadata-binding-audit.json) from
the single main checkout; that record no longer carries the #156 worktree paths.
