# Consolidate the TranSiGen module comparison in CPST Step 05

**Status** Accepted design, 2026-11-28.

Promote Step 06.4 to `cpst/05_tsg_module`. Retire the active Step 05 loss
studies and Step 06.1–06.3, including the width extension. Pin their
reports and replay routes to Git revisions and retain any products still
needed for replay or downstream use. This supersedes the corresponding
rerun requirements in [ADR 0108](0108-cpst-loss-and-shared-evaluation.md).
Steps 07–10 retain their accepted configurations.

Compare `XX-XX`, `XX-XXΔ`, `X-XΔ`, `XΔ-XΔ` and `XΔ-XΔΔ` on the shared
five cold-drug folds and seeds 17, 42, 43, 2026 and 2027. Register 125 new
fits. Keep the Round 4 data, conditioning and architecture definitions.
Use `L = 0.1 K + MSE_X1 + MSE_Δ + MSE_Δpred` for training and validation,
with the full validation loss governing checkpoint selection.

For treated-expression reconstruction, compute reconstructed signed Δ as
reconstructed X2 minus measured X1. With the same X1 in the observed Δ and
the same reduction, the resulting squared error equals the original X2
reconstruction squared error algebraically. Measured-Δ decoders supply
the reconstruction directly. The control-only arm has zero reconstruction
Δ loss and only the control prior KL. The other arms retain their control
and measured-branch prior KLs and predicted-to-measured alignment KL.
Use unit predicted Δ MSE weight throughout.

## Comparison figures

Amend the paired figure rule in [ADR 0109](0109-cpst-seven-metric-comparison-figures.md).
For each metric, compute the arithmetic mean of all valid displayed baseline
fold-and-seed scores. Subtract it from every valid sample, including
baseline samples, without matching coordinates to baseline samples.
Show the baseline group with its centered dots, zero mean and sample SD
when defined. Label the axis as a difference from the baseline mean and
document the sample scope. Keep the existing `metrics-paired` file name.

Make this the default for new figures and provide an explicit historical
replay route. Regenerate only the published Step 04 and new Step 05
`metrics-paired` figures, with their CSV data, receipts and verification.
Other published figures replay from their existing CSVs. Independent
same-fold, same-seed contrast statistics keep their definitions.
