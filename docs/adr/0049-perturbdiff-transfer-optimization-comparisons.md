# Compare Step 06 Transfer optimizations at their own loss minima

**Status: accepted (2026-11-09).**

Extend [ADR 0048](0048-perturbdiff-transfer-diagnosis.md) with five paired
Transfer comparisons. Every intervention uses original Transfer as its
reference and changes one optimization at a time. The dose representation
arms are `random` and `log`. The target-stage training arms are `all-half`,
`main-half`, and `adaLN-only`. Each arm is evaluated at the EMA checkpoint
with the lowest author's validation diffusion loss. Training steps may
differ between arms. This supersedes matched-step evaluation for the
extension. The completed ADR 0048 results remain historical evidence.

## Interventions

| Arm | Change relative to original Transfer |
| --- | --- |
| `random` | Restore the three native, same-seed tenv positive-dose embeddings after loading the original Tahoe bridge. |
| `log` | Use one continuous log-dose encoder shared across Tahoe and tenv, with a separate 0 µM control branch. Refit Tahoe under the same source data and training recipe, then transfer the loss-selected source EMA into tenv. |
| `all-half` | Set the base learning rate of all trainable tenv parameters to `1e-4`. |
| `main-half` | Set the base learning rate of the Cross-DiT denoiser to `1e-4` and retain `2e-4` for the external conditioning encoder. Keep the full model trainable. |
| `adaLN-only` | Update only the existing treated and control adaLN modulation projections in all DiT blocks and both output heads at `2e-4`. Freeze the remaining parameters. |

The three target-stage training arms keep the original dose mapping.
Evaluate all five interventions on Step 05 `cold-cm-drug` fold 0 at seeds
17, 42, and 2026, for 15 target fits. All cells run regardless of the
former screening threshold. Reuse the corresponding original Tahoe and
Transfer checkpoints for each seed. The three `log` source refits are part
of that intervention. The earlier `random` and `all-half` fits used the
`cold-drug` Tahoe source and remain evidence for that regime.

Use the frozen Step 05 `cold-cm-drug` memberships and source-exposure
policy. Tenv holds out drug parent-connectivity blocks, and Tahoe fitting
can include those blocks. See [ADR 0047](0047-perturbdiff-tahoe-transfer-to-tenv.md).

## Evaluation and figure

Use the [Stage 2 metric contract](0042-tenv-model-benchmark-redesign.md#evaluation-and-computation).
For each dose, macro-average AP@100 over the drug–dose groups within each
direction, then give the two directional means equal weight. Label this
metric `AUPRC (mean AP@100)`. PCC is Pearson correlation across all genes
on each group's predicted and true Δ vectors, macro-averaged by dose.

Score the same frozen tenv validation groups in every arm. Publish one
final 2 × 2 figure. Columns are AUPRC and PCC. Rows are dose representation
and learning rate. Each panel plots the metric by tenv dose, with separate
colors for original Transfer and the interventions assigned to that row.
Average the three seed-level values equally at each dose and retain the
individual values for replay. The three seeds are repeat fits on one split.
