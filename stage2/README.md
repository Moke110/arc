# Stage 2 perturbation prediction

## Selected gene axis

**Denoise-DEG 5k (`denoise-5k`)** is the Stage 2 prediction axis.
The [Step 04 matched-control comparison](04_tsg_axis/README.md#gene-axis-selection)
supports this choice based on held-out perturbation Delta prediction.
Across six axes and three cold-drug folds, Denoise-DEG 5k had the highest
mean Δ PCC (0.4066) and the lowest mean Δ RMSE (0.2596) and Δ MMD (0.1793).

Use the [ordered 5,000-gene list](02_gene_axis/02.2_denoise-deg/output/release/gene_axis_5000.csv)
and its [Step 03.7 tenv training release](03_data_train/03.7_tenv_denoise-5k/README.md).

## TranSiGen Delta loss comparison

[Step 05](05_tsg_loss/README.md) studies the Delta loss form and weight on the
selected axis with the same matched plate controls and frozen cold-drug folds.
[ADR 0088](../docs/adr/0088-transigen-delta-loss-comparison.md) defines the
objective grid and its native anchor.

## Local history

[ADR 0087](../docs/adr/0087-stage2-active-steps-and-history.md) defines the archive
cutover. The materialized checkout keeps superseded raw and training releases,
prior comparisons, the archived Step 05–09 studies and rollback materials in
`stage2/history`. The local `archive-manifest.json` records original
paths, archived paths and replay revisions. Normal readers, builds and
comparisons use the active releases. In an archive-isolated checkout with the
active inputs materialized, run `python stage2/04_tsg_axis/execution/verify_active_inputs.py`
to verify the Step 01–03 readers, projections, and fold population.

The [cutover smoke receipt](04_tsg_axis/log/323-archive-isolated-smoke.json)
records the archive-isolated acceptance.

The [pinned archive cutover instructions](https://github.com/Moke110/arc/blob/99078e6f/stage2/README.md#deploy-the-adr-0087-cutover)
record the migration procedure for its original deployment revision.
