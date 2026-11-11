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
| `output/historical-training-release/` | historical v1 truth-fit release (includes the rejected Tahoe payload) |
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

Results reused from the accepted records (this ticket re-verified digests and
path bindings rather than repeating the scientific computation, except where a
bounded gate is named in [`log/156-gates.md`](log/156-gates.md)):

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
  input: its source replay misattributed and dropped counts. Its payload is kept
  as forensic evidence and is never read.
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
The sealed v1 payload is retained as forensic evidence; it is not repaired,
relabeled or claimed equivalent to the current release.

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
