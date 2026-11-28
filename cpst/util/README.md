# CarPS-Transcriptome shared utilities

`cpst/util` provides the identities, frozen features, training and scoring code,
and figure renderers used across CarPS-Transcriptome. Steps 01–03 own the raw
counts, ordered prediction axes and training releases. Steps 04–10 own their
study populations, calibration, score aggregation and reports.
[CarPS-Transcriptome steps](../README.md) · [Utility execution guide](AGENTS.md)

## Shared input and release APIs

Add `cpst/util/output/api` to the Python import path to import `cpst_util`.

| Module | Entry points and contracts |
| --- | --- |
| [`identity`](output/api/cpst_util/identity.py) | `load_gene_table(directory=None)` verifies the indexed table and reads its 18,164 genes in frozen row order. `GeneIdentity` resolves global and source-specific identifiers and rejects ambiguous mappings. [Gene-identity manifest](output/gene-identity/manifest.json) |
| [`representations`](output/api/cpst_util/representations.py) | `RepresentationTable(source="tenv")` or `RepresentationTable(source="tahoe_bulk")` checks the source-specific manifest and its files. `ordered_parents(observations)` collapses repeated records to unique desalted parents in InChIKey order. `table.verified_parent_keys(parents)` checks structures, and `table.vectors(arm, keys)` returns vectors in the caller’s key order, including duplicates. |
| [`ecfp4`](output/api/cpst_util/ecfp4.py) | `open_ecfp4(directory, source_manifest=...)` verifies the 2,048-bit binary fingerprints and parent order, with an optional original manifest binding. `publish_ecfp4(input_directory, candidate_directory)` copies the verified fingerprint and parent table into a separate candidate and records their provenance. |
| [`seeds`](output/api/cpst_util/seeds.py) | `read_seeds(count=None)` verifies the [seed table](output/tables/seeds.csv) and reads `(17, 42, 2026)` in order. Supply `count` to take a prefix. |
| [`tenv_folds`](output/api/cpst_util/tenv_folds.py) | `open_tenv_folds()` verifies the shared release for Steps 04–11. Select a protocol by its versioned ID and align each fold to the ordered treatment table with `fold.positions_for_treatments(treatments)`. |
| [`metrics`](output/api/cpst_util/metrics.py) | `vector_metrics(truth, prediction)` computes PCC, Spearman, RMSE and R² on paired finite vectors. Constant truth leaves both correlations and R² undefined. Constant prediction with varying truth gives zero for both correlations. `signed_jaccard(truth, predicted, gene_indices)` returns combined and direction-specific Jaccard over the top 100 genes per direction. |
| [`publication`](output/api/cpst_util/publication.py) | `publication_lock(destination)` coordinates writers. Hold it while `publish_candidate` builds and verifies a candidate and `replace_release` installs it. |
| [`derivation`](output/api/cpst_util/derivation) | `bulk_deg`, `disease_deg` and `deg_rank` share DEG and ranking operations. `r_env` locates the R worker environment and `release` handles derivation manifests, verification and publication. |

The gene reader defaults to `cpst/util/output/gene-identity/`. The representation
reader defaults to the locally materialized `cpst/util/output/representations/`
for tenv and its `tahoe/` subdirectory for Tahoe. Either reader accepts an
explicit directory. The representation reader uses its built-in source manifest
pin unless `manifest_sha256` is supplied. The manifest pins the input generation,
parent identities, conformer ledger and seven sets of frozen vectors.
[Representation readers and manifest pins](output/api/cpst_util/representations.py) · [Step 08 release evidence](../08_cpst_chem-emb/log/364-formal.json)

## Checkpoint disposition after acceptance

The `cpst_util.checkpoint_disposition` API is for a numbered study with a
completed full acceptance, a pinned acceptance receipt and an independently
verified list of this rerun's eligible run manifests. Keep every weight
required by downstream readers or remaining exact replay out of that list.
Use the numbered study's full verifier as `verify_full`, its numerical and
figure verifier as `verify_retained`, and its accepted
manifest digests as `runs`. `dispose_checkpoints` checks those original
manifests, hashes every declared output of the selected runs and removes only
their `checkpoint.pt` files after the retained verification passes. It stages
weights on the same filesystem during verification and records file sizes and
digests in a versioned receipt without changing the manifests. Preserve an independently trusted
pin for the completed disposition receipt.

`verify_prediction_archive` checks all remaining declared output bytes for the
disposed runs against their original manifests and the pinned completed
receipt. A numbered archive verifier must also verify unaffected runs and
independently rescore predictions and replay figures. The existing full
verifiers still require checkpoints. If disposal stops with a pending receipt,
archive verification fails. Preserve the pending record for investigation.
[Checkpoint disposition API](output/api/cpst_util/checkpoint_disposition.py)

## Shared tenv evaluation

Use `open_tenv_folds().protocol(COLD_DRUG)` for Steps 04–09 and Step 11.
The release fixes five parent groups and seeds `17, 42, 43, 2026, 2027`.
Use `open_tenv_folds().protocol(COLD_DOSE_DRUG)` for Step 10. Its nine folds
cross three historical parent groups with 0.2, 1 and 10 µM in parent-group-major
order. Each fold exposes the held group and dose and both excluded strips.
The ordered training seed roster is `17, 42, 2026`.

Bind new runs to `protocol.name`, `protocol.release_sha256`, `fold.number`,
the training seed and the axis release identity. Use `protocol.coordinates`
to register the full ordered fold-and-seed roster. The positions reader maps
the frozen plate, parent and dose identities to the caller’s matrix row order
and requires the complete population. The returned keys include `train`,
`valid` (validation), and `test`. The [manifest](output/tenv-folds/v1/manifest.json)
records the split seed, roster, ledger and partition digests, and eight
source-axis digests.
Historical Step 03 folds are available through the original readers and
[pinned Step 03 release](https://github.com/Moke110/carps/tree/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/03_data_train).

To reproduce the ledger from all eight materialized tenv axes, run from the
checkout root in `arc` with a fresh candidate directory.

```bash
PYTHONPATH=cpst/util/output/api python cpst/util/execution/build_tenv_folds.py \
  --training cpst/03_data_train --output tmp/tenv-folds-replay
diff -r cpst/util/output/tenv-folds/v1 tmp/tenv-folds-replay
```

## TranSiGen training and evaluation

[`training/transigen.py`](training/transigen.py) owns registration, epoch
selection, minimum-validation checkpoints, prediction and exact replay. It also
provides the architecture KL diagnostics. New runs record the architecture,
objective and coefficients, frozen chemical representation, fitted dose transform
and conditioning, ordered axis digest, split, fold and seed. Numbered run commands
bind the study configuration and execute each phase in a separate process.

| Studies | Shared runtime and accepted evidence |
| --- | --- |
| 04–06 | [`eval/fixed_prediction.py`](eval/fixed_prediction.py) verifies held-out condition and gene order, paired predictions and signed Δ arrays before scoring. [`eval/shape.py`](eval/shape.py) pairs architecture scores with the registered control on the same fold. The consolidation checked 285 installed registrations against the original Git blobs and matched 343 smoke and formal configurations between the original and shared entry points. Representative CUDA fits matched the original checkpoints, predictions and metrics. [TranSiGen verification](training/verification.md) · [Formal receipt checksums](training/formal-receipts.sha256) |
| 07 | [`training/cpst_model.py`](training/cpst_model.py) exposes `CPSTAdapter` for the fixed `XΔ-XΔ` architecture used by current studies. The KPGT default preserves the Step 06.4 initialization and forward path. The alpha release verified its source pins and nine Step 06.4 checkpoint replays. [Alpha release](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/07_cpst_alpha/README.md) · [Pinned alpha release manifest](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/07_cpst_alpha/output/release.json) |
| 08–10 | [`training/chemical.py`](training/chemical.py) aligns frozen chemical vectors to pooled conditions. [`training/registered.py`](training/registered.py) locks runs and validates their phase receipts. [`eval/axis_delta.py`](eval/axis_delta.py) scores the full-axis comparisons in Steps 09 and 10. The acceptance checked the original source and scientific configuration of 207 fresh registrations across Steps 08–10, with reused fits bound to predecessor receipts. [Shared runtime verification](training/cpst-verification.md) · [Receipt checksums](training/cpst-receipts.sha256) |

Step 08’s formal rule retained KPGT as the working default. Changing it required
one representation to lead both primary metrics and each mean paired improvement
to exceed its sample SD. UniMol led both metrics, and its mean Δ PCC improvement
fell below its sample SD. ADR 0100 subsequently adopted frozen UniMol for
CarPS-Transcriptome studies starting with Step 09. [Chemical representation results](../08_cpst_chem-emb/README.md#formal-results) · [Successor decision](../../docs/adr/0100-carps-vae-unimol.md)

Step 09 selected Denoise-DEG 5k on the full prediction axis. Step 10 selected
`no-dose` on the cold-dose-drug intersection. Step 09’s historical common-axis
scores and Step 10’s fold calibration belong to their own reports.
[Prediction axes](../09_cpst_axis/README.md#formal-results) · [Dose inputs](../10_cpst_dose/README.md#formal-results)

Current registrations and receipts use `cpst` paths. The release migration maps
reissued registration digests to the originals embedded in unchanged checkpoints.
The [runtime lookup](output/api/cpst_util/registration-hashes.json) holds that
projection and the accepted Denoise-DEG manifest reissue. The provenance reader
checks historical source digests at each experiment’s pinned Git revision. [Migration inventory](log/385-release-migration.json) · [Pinned Step 06.4 replay](../06_tsg_module/06.4_round4/README.md#replay-a-sealed-step-064-fit)

## Comparison figures

[`comparison`](output/api/cpst_util/comparison.py) renders `learning-curve.png`,
`metrics.png` and `metrics-paired.png` from the published CSVs. Its paired scores
subtract the reference within the same fold and seed. Studies compute
JAC-top200 per held-out condition before aggregating at their established
boundaries. Both metric figures use the seven-panel layout in ADR 0109.
[Current panel layout](../../docs/adr/0109-cpst-seven-metric-comparison-figures.md)
The ADR 0109 migration verifies the prior CSV digests, retains the scored
rows for the current metrics, and renders both layouts from those rows. The
receipt names the predecessor Git revision and records its receipt digest.
The Step 09 and 10 receipts also bind the projected scientific tables.
[Original comparison design](../../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md)

| CSV under `figure/data/` | Columns |
| --- | --- |
| `learning-curve.csv` | `arm,fold,seed,epoch,train,validation` |
| `metrics.csv` | `arm,fold,seed,metric,value` |
| `metrics-paired.csv` | `arm,fold,seed,metric,reference,arm_value,reference_value,value` |

Use `comparison.METRIC_NAMES` for metric IDs and leave undefined values blank.
`write_paired_metrics(data, reference)` derives the paired CSV.
`render_comparison(data, output, arms=..., reference=..., title=...)` takes an
ordered mapping from arm IDs to `(label, hue, depth)`. The renderer checks fit
coverage, epoch order and pair values. `render_receipt(stage, step, output,
renderer, dependencies=...)` runs the renderer twice and records the input CSV,
PNG and renderer hashes. The shared figure gate covers 45 PNGs across the 13
comparisons, Step 02 and the palette. [Figure receipt API](output/api/cpst_util/comparison_publication.py) · [Figure verifier](figure/script/verify_report_figures.py)

## Plotting

The palette has nine hues and four depths. `palette.color_pair(hue, depth=2)`
returns the main and secondary RGBA colors. Hue is in degrees and API depths run
from 0 (darkest) to 3 (lightest). The preview labels the depths 1 to 4. Encode
the first grouping variable with hue and the second with depth. Replicates
inherit their group’s color. [Palette preview](figure/palette.png) · [Plotting design](../../docs/adr/0093-stage2-shared-plotting-and-in-place-figures.md)

[`plotting`](output/api/cpst_util/plotting.py) provides CSV reading, figure saving
and Matplotlib drawing functions. Study renderers control panel layout and labels.

- `bar` accepts a single hue and depth for the series or one of each per bar.
  `grouped_bar(axis, x, groups, width=.8)` centers the supplied series
  `(heights, hue, depth, label)` over each x position.
- `replicate_bar` shows the mean, one dot per defined replicate and sample SD
  when at least two values are defined. `replicate_errorbar` provides the same
  dots and errors on non-bar charts. `replicate_dots` offsets replicates
  symmetrically around the group center.
- `learning_curve(axis, runs)` draws one train and one validation line from
  each run’s `(epochs, train, validation)` values. Call it per model subplot.
- `focus_y(axis, groups, reference=...)` covers finite observations, sample SD
  and the supplied reference. Use the same groups in panels sharing a scale.
  Comparison metric panels pass zero as the reference.
- `line`, `scatter`, `errorbar`, `range_band`, `bubble` and `histogram` use one
  hue and depth per series. `heatmap` uses four depths of one hue and a numeric
  color scale. Pass existing bin edges to `histogram` and `xerr` or `yerr` to
  `errorbar`.

Run `python cpst/util/figure/script/render_palette.py` from the repository root
to redraw the preview. Supply `--output PATH` to write a separate PNG.
[Palette renderer](figure/script/render_palette.py)

## Verify the utility APIs

The shared API tests cover tracked metadata and small fixtures. The utility
execution guide gives the commands for these tests and the study-specific
training and scoring test groups. [Utility test commands](AGENTS.md#verify-the-utility-apis)

Use the materialized release checkout for full receipt and scientific verification.
[Steps 04–06 verification](training/verification.md) · [Steps 07–10 verification](training/cpst-verification.md)
