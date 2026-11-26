# Compare TranSiGen prediction readouts, conditioning heads, and latent residuals

**Status:** accepted (2026-11-26).

Use `06_tsg_shape/06.2_round2` to compare prediction-path changes on the
Step 06.1 original architecture. Hold the frozen tenv cold-drug folds,
seed 42, Denoise-DEG 5k, measured same-plate DMSO, frozen KPGT, standardized
numeric dose, and training settings fixed. Reuse the three sealed original
fits from Step 06.1 after verifying their conditions, configuration, and
evaluation provenance.

| Coordinate | Change from `original` |
| --- | --- |
| `original` | Accepted Step 06.1 original architecture and Step 05.2 objective. |
| `output-delta` | Add an independent Delta decoder with the original decoder's hidden structure and a signed linear output. Feed it the sample from the predicted treated latent distribution. Keep both expression reconstructions and the Step 05.2 objective. |
| `output-delta-nox2` | Use the original treated decoder to decode the predicted treated latent sample into signed Delta. Remove treated expression reconstruction. Keep the measured treated encoder for the prediction KL and keep all three KL terms. Replace the treated reconstruction MSE and prediction Delta MSE with twice the prediction Delta MSE. |
| `z-interact_MLP` | Concatenate the 100 latent, 400 drug, and 32 dose features. Use a 532→145→200 ReLU MLP and split its outputs into predicted mean and log variance. |
| `z-interact_FiLM` | Use a 432→136→200 ReLU MLP on the drug and dose features to generate 100 scale offsets and 100 shifts. Apply `(1 + gamma) * Z1 + beta` and use a 100→200 linear projection split into predicted mean and log variance. Zero-initialize the modulation MLP's final layer for identity modulation at initialization. |
| `shared-z` | Share the full encoder trunk and its BatchNorm parameters and running statistics across control and treated expression. Keep separate Gaussian parameter projections. |
| `z-res` | Add the sampled control latent vector Z1 to the original linear prediction head's mean. Keep its log-variance projection. |
| `shared-z-res` | Apply the same shared encoder trunk and prediction-mean residual together. |

For both direct Delta arms, construct predicted treated expression as
`X2_hat = X1 + Delta_hat` using the measured matched control. The signed
linear output keeps both directions of Delta representable. In
`output-delta-nox2`, the prediction-space expression MSE
`MSE(X2_hat, X2)` equals `MSE(Delta_hat, X2 - X1)` algebraically.
The prior treatment autoencoder reconstruction term is dropped. Use
`0.1 * K + MSE_X1 + 2 * MSE_Delta + 1.75 * L1_Delta` for training and
checkpoint selection in that arm, with all three native KL terms in `K`.
Treat comparisons involving either direct Delta arm as effects of the
registered architecture and objective changes.

Preserve the 32 repeated standardized dose features in both interaction
arms. The MLP and FiLM prediction heads have 106,485 and 106,488 trainable
parameters, respectively, compared with 106,600 in the original head.
For the shared trunk, preserve the normal control-then-treated training
calls, with a single set of BatchNorm statistics. An identity skip on a
linear mean head is absorbable into its weights. Its effect in this
comparison reflects parameterization, initialization, and regularization.

Train and independently replay the 21 new fits. Use held-out Delta PCC and
Delta RMSE as primary outcomes and publish all nine metrics with per-fold
scores, equal-fold means, and within-fold differences from `original`.
For each score or difference panel, draw all eight groups on the same axes
in the order above using the shared plotting API. Publish one learning-curve
panel per group and the CSVs needed to reproduce every figure. Identify the
objective difference in the `output-delta-nox2` learning curves.

[Matched controls](0086-tenv-axis-matched-plate-control.md) ·
[Step 05.2 objective](0090-transigen-mse-l1-sweep.md) ·
[Step 06.1 baseline](0094-transigen-shape-round1.md) ·
[Sealed Step 06.1 fits](0095-transigen-shape-encoder-width-extension.md)
