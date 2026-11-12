# Stage 2 shared foundations

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
