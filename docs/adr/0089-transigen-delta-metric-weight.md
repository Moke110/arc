# Vary the TranSiGen Delta loss coefficient independently

**Status:** accepted (2026-11-25).

Place the loss comparison in `05_tsg_loss/05.1_loss_metrics`. Hold both
reconstruction MSE coefficients at 1 and the sum of the three native KL
terms at 0.1. Apply coefficients 0.75, 1 and 1.25 directly to the Delta
loss to isolate its weight from reconstruction weighting. Cross those
coefficients with MSE, Cos, PCC and L1. MSE at coefficient 1 calls the
native objective and provides the baseline for score differences.

Keep the matched control, fixed `denoise-5k` inputs, three frozen cold-drug
folds, seed 42, loss reductions and held-out score definitions from
[ADR 0088](0088-transigen-delta-loss-comparison.md). Replace the prior
Step 05 study results with the verified Step 05.1 matrix and the raw score,
score difference and learning-curve figures. This replacement has no
separate archived result set.
