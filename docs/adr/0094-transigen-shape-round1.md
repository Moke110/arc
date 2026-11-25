# Compare TranSiGen expression widths and depths in Stage 2 Step 06.1

**Status:** accepted (2026-11-25).

[ADR 0095](0095-transigen-shape-encoder-width-extension.md) extends this
registered eight-coordinate comparison.

Use `06_tsg_shape/06.1_round1` to compare expression encoder width,
expression latent dimension, and expression encoder and decoder depth.
Changing each independently isolates it from the joint width scaling studied
in ADR 0081. The combined depth coordinate also completes the encoder-depth
by decoder-depth comparison.

| Coordinate | Encoder hidden widths | Decoder hidden widths | Expression latent dimension |
| --- | --- | --- | ---: |
| `widen_z-100` | 1200 | 800 | 100 |
| `widen_z-110` | 1200 | 800 | 110 |
| `widen_z-120` | 1200 | 800 | 120 |
| `widen_z-130` | 1200 | 800 | 130 |
| `widen_enc` | 6000 | 800 | 100 |
| `deepen_enc` | 1200, 1200 | 800 | 100 |
| `deepen_dec` | 1200 | 800, 800 | 100 |
| `deepen-both` | 1200, 1200 | 800, 800 | 100 |

Apply expression changes to the control and treated branches, preserving
their independent parameters and the sharing of the treated decoder across
reconstruction and prediction. Added blocks use the native linear,
BatchNorm, ReLU and dropout sequence. Apply each latent dimension jointly to
control, measured treated and predicted treated latent vectors and their
distribution heads. The conditioning head input is the latent dimension
plus 400 drug features and 32 dose features.

Fit all eight coordinates on each of the three frozen tenv cold-drug folds
with seed 42, Denoise-DEG 5k and measured same-plate DMSO throughout.
Use Step 05.2's objective at Delta L1 weight 1.75, unit expression and
Delta MSE weights, and weight 0.1 on the three summed native KL terms.
Retain Step 05.2's inputs, optimizer, sampling seeds and minimum validation
objective checkpoint rule with patience 30 and a 1,000-epoch cap. The KL
sum covers more latent dimensions as Z grows. Log the total and the
per-latent-dimension KL for interpretation. Refit the original architecture
at `widen_z-100` and check it against the accepted Step 05.2 weight-1.75 runs.

Plot all nine held-out Delta metrics with fold scores, equal-fold means
and paired differences from the fresh control. Use metric-specific
direction and preserve undefined-score handling from Step 05.2. Report
results descriptively without ranking configurations, choosing a winner or
applying an overall improvement rule. Reuse the Stage 2 plotting API and
publish the inputs for CSV-only figure replay.

- [ADR 0062](0062-stage2-accepted-products-and-reproducible-figures.md)
- [ADR 0081](0081-transigen-width-and-dose-comparison.md)
- [ADR 0090](0090-transigen-mse-l1-sweep.md)
- [ADR 0093](0093-stage2-shared-plotting-and-in-place-figures.md)
