# Compare Step 06 Transfer optimizations at their own loss minima

**Status: proposed (2026-11-09).**

Extend [ADR 0048](0048-perturbdiff-transfer-diagnosis.md) with five paired
Transfer comparisons. Every intervention uses original Transfer as its
reference and changes one optimization at a time. The dose representation
arms are `random` and `log`. The target-stage training arms are `all-half`,
`main-half`, and `adaLN-only`. Each arm is evaluated at the EMA checkpoint
with the lowest author's validation diffusion loss. Training steps may
differ between arms. This supersedes matched-step evaluation for the
extension. The completed ADR 0048 results remain historical evidence.

Publish one final 2 × 2 figure. Columns are AUPRC and PCC. Rows are dose
representation and learning rate. Each panel plots the metric by tenv
dose, with separate colors for original Transfer and the interventions
assigned to that row.

Before the spec, resolve the evaluation population and confirmation gate,
the AUPRC definition, the shared log-dose encoder's fit scope, and the
parameter groups for `main-half` and `adaLN-only`.
