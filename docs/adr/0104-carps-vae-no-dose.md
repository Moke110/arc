# Select no-dose for subsequent CarPS-VAE work

**Status** Accepted, 2026-11-27.

Use the Step 10 `no-dose` configuration for subsequent CarPS-VAE work.
The prediction head takes the control latent and the projected frozen UniMol
representation. Measured dose still identifies treatment conditions and
defines the held-out partitions in the Step 10 evaluation.

Across three cold-dose-drug folds and seeds 17, 42 and 2026, `no-dose` had
the highest mean full-axis held-out Δ PCC (0.3913 ± 0.0235) and the lowest
mean Δ RMSE (0.2472 ± 0.0046). The best dose-input arm, `Zdose-MLP`, scored
0.3700 ± 0.0324 and 0.2511 ± 0.0041. Its paired differences relative to
`no-dose` were −0.0212 ± 0.0254 in Δ PCC and +0.0039 ± 0.0056 in Δ RMSE.
Values are means ± sample SD across nine fold and seed coordinates.
[Step 10 formal results](../../stage2_CarPS-VAE/10_cpsv_dose/README.md#formal-results)

The evidence covers simultaneous generalization to unseen drug parents and
doses on the frozen tenv Denoise-DEG 5k release. Each held dose is paired
with a different group of held parents, so foldwise differences combine
both effects. [ADR 0103](0103-carps-vae-dose-inputs.md)
