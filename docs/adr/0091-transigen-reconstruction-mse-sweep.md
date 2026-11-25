# Vary the shared TranSiGen reconstruction MSE weight

**Status:** accepted (2026-11-25).

Use `05_tsg_loss/05.3_mse-l1_mse-x` to fit
`0.1 * K + MSE_Delta + 1.75 * L1_Delta + alpha * (MSE_X1 + MSE_X2)`.
`K` is the sum of the three native KL terms. Sweep the shared reconstruction
MSE weight over 0, 0.5, 0.75, 0.9, 1, 1.1, 1.25, 1.5 and 2 on the three
frozen cold-drug folds with the protocol from
[Step 05.2](0090-transigen-mse-l1-sweep.md).

Use alpha 1 as the signed held-out score-difference reference so the study
measures changes from the established Step 05.2 objective at Delta L1 weight
1.75. Refit the reference on all three folds and preserve the predecessor's
objective arithmetic at this coordinate to verify exact agreement.
Alpha 0 removes both reconstruction MSE terms while retaining both Delta
terms and all three KL terms.
