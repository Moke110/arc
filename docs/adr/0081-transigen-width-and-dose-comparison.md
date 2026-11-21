# Compare TranSiGen width and dose encoding on tenv

**Status** — accepted (2026-11-21).

Stage 2 Step 05 crosses no dose, latent numeric, latent embedding, and FiLM
with TranSiGen width multipliers 1, 1.5, 2, and 2.5. The same multiplier
scales the internal expression encoder and decoder, latent, chemical
embedding, and applicable dose-conditioning vector dimensions. The frozen
KPGT input and 5,000-gene expression interfaces stay fixed. One rule at
every width isolates capacity from dose-encoding changes.

Reuse the [ADR 0080](0080-tenv-kpgt-dose-round1.md) tenv freeze, matched
controls, and three cold-drug-dose folds. Train all 48 coordinates from
scratch using Step 04's seed and training recipe, selecting each checkpoint
by minimum treated-expression validation MSE. The twelve original-width
coordinates can be checked against their Step 04 counterparts. Move the
existing transfer step to Step 06. Its existing protocol is independent
of the Step 05 comparison.

Compute Δ Pearson, Δ Spearman, Δ RMSE, and Δ R² with equal sample weight
within each fold. Summarize each width and encoding with an equal-fold
mean. The bar chart shows individual fold scores and their sample
standard deviation. In the overall balloon chart, color and text show
the mean and diameter increases with its z-score relative to all 48 run
scores for that metric. The metrics keep their original direction.
The per-dose chart uses the same sample-equal score as the corresponding
fold dot and shares color and diameter scales across widths within each
metric. Each dose corresponds to one fold, so dose and fold effects
cannot be disentangled in these panels.
