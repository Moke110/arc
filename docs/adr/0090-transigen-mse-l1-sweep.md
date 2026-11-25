# Add a Delta L1 penalty to the native TranSiGen objective

**Status:** accepted (2026-11-25).

Use `05_tsg_loss/05.2_mse_l1` to measure the added penalty's effect while
preserving the native coefficients. Add `alpha * L1_Delta` to the native
TranSiGen objective with unit MSE weights and KL weight 0.1. Use the native
objective's arithmetic and the [Step 05.1](0089-transigen-delta-metric-weight.md)
L1 batch sum across conditions and genes. Both Delta vectors use the measured
same-plate DMSO control. Sweep alpha from 0 through 2 in steps of 0.25 across
the three frozen tenv cold-drug folds on denoise-5k.

Refit all three zero-weight controls using the native loss directly and verify
them against the accepted Step 05.1 native MSE runs. Refitting records the
unweighted Delta L1 alongside the other loss components at every epoch and
gives every coefficient the same audit coverage. Use alpha zero as the
reference for signed held-out score differences. Report the complete sweep
per metric and per fold.
