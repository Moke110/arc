# Sweep bounded TranSiGen Delta similarity loss weights

**Status:** accepted (2026-11-25).

Use `05_tsg_loss/05.4_mse-l1_pcc` for Pearson correlation and
`05_tsg_loss/05.5_mse-l1_cos` for cosine similarity. Fit
`0.1 * K + MSE_X1 + MSE_X2 + MSE_Delta + 1.75 * L1_Delta + alpha * Lmetric`.
Sweep `alpha` over 0, 0.1, 0.2, 0.3, 0.4 and 0.5 in each step.

Compute `Lmetric` as the batch mean of `(1 - clamp(similarity, -1, 1)) / 2`
over conditions, with equal weight per condition in each batch. Compute each
similarity across the ordered gene axis from predicted and measured Delta
against the same measured plate control. Reuse Step 05.1's centering for
PCC, cosine normalization with epsilon `1e-8`, and similarity zero for a
zero-norm vector. Identical, orthogonal and opposite nondegenerate vectors
give losses of zero, one half and one respectively. The metric contribution
is at most `alpha` per batch. Retain the Step 05.3 batch sums for the base
objective and its epoch reporting and checkpoint-selection aggregation of
batch objectives divided by the partition condition count.

Inherit the Step 05.3 model, inputs, folds, random seed, runtime, stopping
settings and reporting protocol. Refit all three alpha-zero folds separately
in each study, using the exact Step 05.3 unit-weight objective for the loss
and gradient. Require exact agreement with Step 05.3 on shared training and
validation traces, selected and stopping epochs, validation and test
predictions, Delta arrays, scores and eligibility. Each study subtracts its
own fresh alpha-zero scores for its signed score differences.

## See

- [ADR 0089](0089-transigen-delta-metric-weight.md)
- [ADR 0091](0091-transigen-reconstruction-mse-sweep.md)
