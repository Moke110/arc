# Stage 2 shared foundations

**Current expression release: ADR 0077, 2026-11-18 (#243–#245).**
Only paired tenv-5k and tahoe_bulk-5k expression is live. Each profile uses its
unchanged G_5000 integer-count total for CP10K and natural log1p. The exact four
ADR0073 gene axes and all disease/fusion/HVG evidence remain unchanged.

[Acceptance](log/acceptance.md) and [the current receipt](log/axis_local_acceptance.json)
pin the paired independent replay, selected totals, preservation, document and
health checks. Earlier measured receipts remain sealed historical full-axis
releases. Stage 2 02/03 frozen inputs, weights and reported metrics retain their
original full-axis scale; this release does not retrain or rescore them.

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

### Complete disease evidence (ADR 0073)

The five category directories under `output/disease-deg/` contain 41 dataset
records: nHF 8, ICM 5, AF 5, genetic_cardiomyopathy 11 and TOX 12. These are
**dataset records, not 41 independent physical studies**: shared control groups,
disease arms and per-donor pooling are explicit in each catalog. All qualifying
declared tiers are processed without a cumulative DEG-count stop. Dated
`source_review` records retain finite primary-source decisions, exclusions,
uncertainty and reentry conditions. See [the source overview](../../disease_data/README.md).

Each source keeps its compact `*.deg.csv` and frozen-axis `*.rank.csv`. The
pass-any-contrast rule is `padj < 0.05` AND `abs(log2FC) >= 2`; non-DEG/missing
source genes have rank `n_source_DEG + 1`, normalized by that source's DEG count.
A zero-DEG/non-informative source is retained with its reason and contributes no
rank term. Evidence thresholds, sample identities and CM/model admission remain
unchanged. Human atrial/ventricular bulk is tier 3, not CM-specific tier 1.

Catalog schema is `stage2-disease-deg-catalog-v2` with the unchanged
`stage2-disease-deg-v2` source rank contract. Catalogs pin source/derived inputs,
contrast/sample units, per-source tier, informativeness, `n_deg`,
`n_deg_in_deg1000`, source review and complete cumulative `tiering.tier_union`.
The collection policy is `all_qualifying_declared_sources`; there is no stop
limit, stopped tier or back-filled missing source tier. Source bytes are immutable
ADR 0030 inputs; declared derivatives stay item-local.

```python
import sys; sys.path.insert(0, 'stage2/01_bases/output/api')
from bases.disease_deg import DiseaseDeg, DiseaseFusion
db = DiseaseDeg('stage2/01_bases/output/disease-deg')
table = db.table('nHF', 'gse270788__nicm')
rank = db.rank('nHF', 'gse270788__nicm')
fusion = DiseaseFusion('stage2/01_bases/output/disease-fusion')
full_evidence = fusion.union('nHF')
retained = fusion.deg1000('nHF')
```

### Tier-first fusion and disease eligibility

For every gene in a disease's complete axis-resolved DEG union,
`source_tier_highest` is the minimum numeric tier among admitted informative
sources where it is actually DEG-positive. Mere measurement confers no tier.
`disease_rank_sum = sum_source(rank / n_source_DEG)` stays a numeric, unweighted
score. Order the full union by highest tier ascending, score ascending, then
original `gene_index` ascending, and assign contiguous ranks `1..n_d`.
`*.union.csv` retains the complete ordering and supporting source IDs;
`*.deg1000.csv` is its first `min(1000, n_d)` rows, with no padding.

| disease | dataset records | complete union | DEG1000 | eligible (`n_d >= 100`) |
|---|---:|---:|---:|---|
| nHF | 8 | 1,785 | 1,000 | yes |
| ICM | 5 | 1,671 | 1,000 | yes |
| AF | 5 | 25 | 25 | no, fewer than 100 |
| genetic_cardiomyopathy | 11 | 4,399 | 1,000 | yes |
| TOX | 12 | 9,589 | 1,000 | yes |

AF exclusion is a measured gate result, not a permanent category exception.
Excluded evidence remains readable and its cross-disease rank/coverage is
unavailable with an explicit reason. TOX is one disease: each of its twelve
admitted dataset designs pools independently supported small-molecule drugs
against matched vehicle. Mechanism is annotation, not a separate disease term.
Its source/drug admission and matched-well audits remain in `output/disease-deg/tox/`.

## Layout

```text
output/
  api/                 # the single public implementation
  gene-identity/       # frozen 18,164-gene identity
  raw_counts/          # native/projected integers
  corpus/              # full-axis integer counts, roles, matched controls
  tenv-hvg/            # descriptive ranking + retained train-only rankings
  disease-deg/         # five categories, 41 source-record DEG/rank pairs
  disease-fusion/      # full *.union.csv + tier-first *.deg1000.csv
  gene-axis/           # candidate table, four axes, composition/coverage, manifest
  view-5k/             # paired tenv + tahoe_bulk X/control/delta payloads
  fold_n_seed/         # retained partitions
  palette/             # shared colors
figure/{data,script}/   # complete plotted CSVs and CSV-only renderers
log/                   # acceptance/replay evidence
gene_axis_fill_rule.docx # portrait A4 human reading aid, not a model input
```

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
## Four exact nested gene axes (ADR 0073)

For each `N in {2000,3000,4000,5000}`, take the first `N/2` tenv-HVG genes,
then unused disease-supported candidates (`f1000 >= 1`) ordered by `f1000`
descending, cross-disease `rank_sum` ascending, `gene_index` ascending. Fill a
shortage only with unused HVGs in HVG order. `f1000` counts eligible diseases
whose DEG1000 contains a gene. Cross-disease `rank_sum = sum_d(rank_d / n_d)`
uses each **complete** disease-union size `n_d`; within DEG1000 use its tier-first
disease rank, otherwise `n_d + 1`, including a true DEG ranked 1,001 or later.
Excluded diseases and tenv-HVG contribute no disease term.

The measured four axes have exact sizes and are nested. Each has half HVG
allocation and half disease candidates, with **zero measured fallback** at all
four sizes; fallback remains implemented and tested. The complete 18,164-row
candidate table has `f1000` in `0..4`. At 5k, DEG1000 coverage is genetic/TOX
1,000 each, nHF 869 and ICM 799. No old result was assumed in the rebuild.

`output/gene-axis/` contains `candidate_rank_sum.csv`,
`gene_axis_{2000,3000,4000,5000}.csv`, `axis_composition.csv`,
`source_topn_coverage.csv` and `manifest.json`. Schema
`stage2-gene-axis-release-v3` and algorithm
`tier-first-deg1000-f1000-four-axes-v3` reject mixed old contracts.
The compatibility coverage **filename** is retained; its explicit long columns
`N,n_axis,group,group_label,numerator,denominator,rate,note` carry current
DEG1000 semantics. Eligible disease denominator is its own DEG1000 size;
tenv-HVG denominator is N; excluded AF rows are unavailable, not zero.

## Paired 5k-local expression and bounded readers

`output/view-5k/` contains float32 tenv (4,478 rows: 4,206 treated/272 vehicle)
and tahoe_bulk (66,718 rows: 65,318 treated/1,400 vehicle). Original identities,
order, roles and matched-control groups are unchanged. First select G_5000 raw
integer counts; sum those counts in int64, then compute
`X = log1p(10000 * counts_5k / selected_axis_library_total)`.
Normalize each vehicle independently before averaging its matched group;
`control` is that mean and `delta = X - control`. Non-positive selected totals
fail with source/profile identity. No full-axis diagnostic is a denominator.

The `stage2-view-v3` manifest, payload metadata and 71,196-row
`library_totals.csv` bind this recipe and the original ascending axis digest.
The reader rejects historical full-axis and mixed-generation payloads. The
2k/3k/4k expression products are retired; their named axes with symbols remain.

```python
from bases.view import open_view
with open_view('stage2/01_bases/output/view-5k') as view:
    tenv = view.source('tenv')
    x, control, delta = tenv.X(0, 64), tenv.control(0, 64), tenv.delta(0, 64)
    controls = tenv.matched_control_rows(42)
    genes = tenv.gene_index
```

Use bounded public windows/batches. Future Stage 2 02/03 input preparation uses
`execution/axis_local/prepare.py` and a fresh empty candidate. It records actual
view schema, recipe, normalization, axis and digests. The original prepare code
is retained for historical fit-identity reconstruction; use the preserved frozen
inputs to replay accepted experiments, whose metrics are on the old scale.

## Readers

Import `output/api/bases` through `execution/public_api.py` or add its parent
to `sys.path` from the repository root. Every product read checks its digests.

| product | reader |
|---|---|
| `raw_counts/` | `bases.raw_counts.open_tenv` / `open_tahoe` |
| `corpus/` | `bases.corpus.open_corpus` |
| `view-5k/` | `bases.view.open_view` (bounded X/control/delta, library totals) |
| `disease-deg/` | `bases.disease_deg.DiseaseDeg` (categories/datasets/table/rank/eligibility) |
| `disease-fusion/` | `bases.disease_deg.DiseaseFusion` (record/union/deg1000) |
| `gene-axis/` | `bases.gene_axis.GeneAxis` (axis/candidate_table/composition/coverage/digests) |
| `gene-identity/` | `bases.identity` |

## Rebuild, replay and figures

Run one heavy job at a time under `execution/capped.sh` and `/usr/bin/time -v`.
Drivers stage outside formal output under `execution/checkpoints/`, verify, then
replace; a failed candidate leaves accepted products untouched. Source-item
recipes are declared in each item's `author_log.json`; preserve ignored bytes
before deleting worktrees. `prepare` of the integrated successor snapshots a
complete materialized input root and independently replayed TOX output. Its
publish/verify actions cover registry categories, fusion, four axes and only the
paired 5k view. This scientific rebuild driver is separate from a scoped
normalization-only publication: use `publish_axis_local.py` for the latter.

```bash
python stage2/01_bases/execution/build_disease_deg.py publish --all
python stage2/01_bases/execution/build_disease_deg.py verify --all
ARC_RSCRIPT=/home/chang/miniconda3/envs/biomni_e1/bin/Rscript python stage2/01_bases/execution/build_tox.py publish
ARC_RSCRIPT=/home/chang/miniconda3/envs/biomni_e1/bin/Rscript python stage2/01_bases/execution/build_tox.py verify
python stage2/01_bases/execution/build_disease_fusion.py verify
python stage2/01_bases/execution/build_gene_axis.py --verify
python stage2/01_bases/execution/build_view.py verify
python stage2/01_bases/figure/script/render_gene_axis_rank_sum_histogram.py --csv stage2/01_bases/figure/data/gene_axis_rank_sum_histogram.csv --png stage2/01_bases/figure/gene_axis_rank_sum_histogram.png
python stage2/01_bases/figure/script/render_gene_axis_source_topn_coverage.py --csv stage2/01_bases/figure/data/gene_axis_source_topn_coverage.csv --png stage2/01_bases/figure/gene_axis_source_topn_coverage.png
python stage2/01_bases/execution/make_gene_axis_docx.py --out stage2/01_bases/gene_axis_fill_rule.docx
```

Registry source fits resolve arc R through `bases.r_env`; TOX uses the documented
biomni R/edgeR/DESeq2 environment. Avoid a global R override across both recipes.
Per-view publish/replay may be run as separate capped public calls to bound process
lifetime; the final combined receipt must validate the entire input/output chain.
The source fits and output hashes in `log/source_derivation_replay.json` are real
model runs, not patched implementation digests.

Figure renderers consume only their tracked CSVs plus the palette; no matrix or
scientific recomputation. The histogram has four axis panels. The relocated
portrait A4 DOCX has five disease subtables (all 41 dataset records) plus separately
labelled tenv-HVG evidence. Use the pinned Python/docx versions in
`log/gene_axis_document_verify.json`; byte replay and every-page rendered QA are
both required. The accepted document has 13 inspected pages. Remove its old
live `output/gene-axis/` copy when promoting a coherent generation.

```bash
script/check_paths.sh
script/check_tracked.sh
script/verify_all.sh   # alone; hashes all external manifests
```

[Log routing](log/README.md) identifies current receipts and sealed predecessors.
Unrelated raw-count/corpus/HVG bytes remain unchanged; no model was trained.

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


## Historical boundary

The errata above and sealed acceptance sections retain original measured facts.
ADRs 0066–0072 describe their own releases; ADR 0073 defines retained disease/axis evidence; ADR 0077 and the current receipt
define the paired 5k-local expression release.

## Scoped publication and recovery

`build_view.py publish|verify` stages/replays only paired 5k expression from
retained corpus integers and G_5000. To install an already independently replayed
candidate, run `execution/publish_axis_local.py install --root <materialized-root>
--candidate <paired-candidate> --recovery <fresh-checkpoint> --record <receipt>`.
The driver verifies both sources and every input pin before moving any live
path, restores moved paths on failure, then retains superseded views outside
`output/`. Recovery instructions are in the driver's module documentation.
Recovery is historical full-axis evidence requiring its matching old reader.
No active command silently creates small expression products.
