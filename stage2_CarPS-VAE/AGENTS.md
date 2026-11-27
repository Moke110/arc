# AGENTS.md — Stage 2

Before working in a step, read its `AGENTS.md` and `README.md`, ADRs 0062 and 0087,
and the scientific ADRs it names. Apply the repository execution-environment rules.

## Active steps

- [util](util/AGENTS.md) owns gene and drug identities, representations, ECFP4,
  seed selection, generic vector metrics, the palette and CSV-only plotting.
  It also owns release publication and locking, shared DEG and ranking kernels,
  R worker resolution, and derivation manifest verification. Read the
  [public API inventory](util/README.md) and reuse compatible operations.
  Keep fold validation, metric aggregation and report tables with the study.
- [01_data_raw](01_data_raw/AGENTS.md) publishes source-native integer counts, coverage
  and match status.
- [02_gene_axis](02_gene_axis/AGENTS.md) owns the ordered HVG-DEG and Denoise-DEG
  prediction axes and their selection evidence.
- [03_data_train](03_data_train/AGENTS.md) publishes eight plate-pooled tenv axes
  and the Tahoe well-level HVG-DEG 5k release.
- [04_tsg_axis](04_tsg_axis/AGENTS.md) compares fixed tenv prediction axes with TranSiGen.
- [05_tsg_loss](05_tsg_loss/AGENTS.md) compares TranSiGen Delta losses and reconstruction weights on the fixed denoise-5k axis.
- [06_tsg_shape](06_tsg_shape/AGENTS.md) compares TranSiGen expression widths, depths and prediction architectures on the fixed denoise-5k axis.
- [07_cpsv_alpha](07_cpsv_alpha/AGENTS.md) owns the fixed CarPS-VAE alpha entry point, KPGT default and exact checkpoint replay.
- [08_cpsv_chem-emb](08_cpsv_chem-emb/AGENTS.md) compares seven frozen chemical representations using the fixed alpha.
- [09_cpsv_axis](09_cpsv_axis/AGENTS.md) compares eight CarPS-VAE prediction axes with UniMol.

Use Step 04's `execution/matrix.py` and `execution/evaluate.py` for the axis comparison.
Steps 05.1–05.5 and 06.1–06.4 reuse the matched-control runtime and metric primitives with their own run and report commands.

## CarPS-VAE configuration

Use the [selected CarPS-VAE configuration](README.md#selected-carps-vae-configuration)
for the model release and its successors. Replay each sealed TranSiGen study
with its pinned protocol. Use the [pinned Step 06.4 route](06_tsg_shape/06.4_round4/README.md#replay-a-sealed-step-064-fit)
when a registration binds the prior source layout.

## Local history

Put retired releases, study steps and rollback material in the Git-ignored
`history/` archive. Preserve sealed bytes and record original paths and replay
revisions in its `archive-manifest.json`. Verify the archive and active entry
points with history unavailable before completing a cutover.

## Step contract
- `output/`: minimal accepted products reused by other steps: corpora, gene/HVG tables,
  evaluator truth, public metric/reader code and reusable model checkpoints. Each has
  one stable path and one owner; public implementations live here and execution imports them.
- `figure/*.png`: human-facing report figures.
- `figure/data/{figure}_{panel}.csv`: complete plotted values for each panel, including
  curve coordinates, groups, uncertainty and reference values where applicable.
- `figure/script/`: reproduce PNGs from those CSVs alone, without matrices, checkpoints
  or recomputing scientific metrics. Record provenance when a CSV derives from output.
  Draw marks with the [shared util plotting API](util/README.md#plotting).
  Plot learning curves as lines with one subplot per model configuration.
  Draw one train and one validation trace for every `(fold, seed)`.
  Use bars or grouped bars for other figures unless the request specifies
  another chart type. Use each `(fold, seed)` as a within-group replicate.
  Draw defined replicate values as dots, bar heights as means and error
  bars as sample SD when at least two values are defined.
  Set y-axis limits near the plotted values when this makes between-group
  differences easier to see. The range may omit zero. Keep the plotted dots,
  complete error bars and relevant references visible. Indicate when a bar
  chart omits zero, and use common limits for directly comparable panels.
  Encode each grouping variable on its own hue or depth axis, up to two
  grouping variables. Replace figures at their canonical paths and discard
  superseded render bundles.
  Keep a single current set of figures, plotted data and renderers under
  each `figure/`, using only `data/` and `script/` as subdirectories.
  [Plotting decision](../docs/adr/0093-stage2-shared-plotting-and-in-place-figures.md)
  [Current figure acceptance](README.md#verify-the-current-figures)
- `execution/`: reproducible build/audit code, configurations, tests and disposable
  working intermediates. Keep necessary reconstruction dependencies available.
- `log/`: execution/audit records, environments, seeds, acceptance results, hashes,
  resource measurements, migration mappings and compact failure evidence.
- `AGENTS.md`: agent execution/reproduction workflow; `README.md`: human-readable
  accepted results, limitations and links to figures and reusable products.

## Acceptance and replacement

Build candidates outside formal output, verify structure/data constraints, scientific
criteria, downstream reads and CSV-only figure reproduction, then replace stable paths.
Record metric tolerances for stochastic results and hashes for file integrity. A failed
candidate leaves the accepted products available. Keep accepted final products active and
archive superseded scientific payloads after acceptance and reference migration.
Treat external `source/` trees as read-only.
