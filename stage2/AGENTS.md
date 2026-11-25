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
- [03_data_train](03_data_train/AGENTS.md) publishes six plate-pooled tenv axes
  and the Tahoe well-level HVG-DEG 5k release.
- [04_tsg_axis](04_tsg_axis/AGENTS.md) compares fixed tenv prediction axes with TranSiGen.
- [05_tsg_loss](05_tsg_loss/AGENTS.md) compares TranSiGen Delta losses and reconstruction weights on the fixed denoise-5k axis.
- [06_tsg_shape](06_tsg_shape/AGENTS.md) compares TranSiGen expression widths and depths on the fixed denoise-5k axis and selected objective.

Use Step 04's `execution/matrix.py` and `execution/evaluate.py` for the axis comparison.
Steps 05.1–05.5 and 06.1 reuse the matched-control runtime and metric primitives with their own run and report commands.

## TranSiGen defaults

For new Stage 2 TranSiGen experiments, use the ordered Denoise-DEG 5k
(`denoise-5k`) axis. Use measured same-plate DMSO as model input `X1` in
training, validation and test. Subtract that same control from both predicted
and measured treated expression to compute Δ.

Use `L = 0.1 * K + MSE_X1 + MSE_X2 + MSE_Δ + 1.75 * L1_Δ` with the three
native KL terms summed in `K` and unit weights on both reconstruction MSEs.
Use each pinned study protocol when replaying recorded results.
[Selected configuration and evidence](README.md#selected-transigen-configuration)

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
  Aggregate fold and seed replicates within each experimental setting. Never use
  fold or seed for plot axes, facets, legend entries, or visual groupings.
  Use a named `figure/<successor>/` subtree when a receipt seals the predecessor's
  paths. Bind the successor inventory to its own acceptance receipt.
  [ADR 0093](../docs/adr/0093-stage2-shared-plotting-and-in-place-figures.md)
  authorizes an in-place refresh of the 22 specified report figures. Use the
  new rendering receipts for those paths and the pinned Git revision for
  original PNG replay. [Current figure acceptance](README.md#verify-the-current-figures)
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
archive superseded payloads after acceptance and reference migration.
Treat external `source/` trees as read-only.
