# Step 06.1 TranSiGen expression architecture

Step 06.1 tests how encoder width, latent dimension, and encoder and decoder
depth affect held-out expression prediction. The `original` architecture uses
1,200 encoder units, 800 decoder units, and a 100-dimensional latent. The first
matrix added latent widths 110, 120, and 130, encoder width 6,000, and a second
hidden layer on either side or both. The extension added encoder widths 2,000
and 4,000, bringing the comparison to ten architectures. Both expression
branches use each registered shape with independent parameters. The treated
decoder is shared between reconstruction and prediction.
[Initial design](../../../docs/adr/0094-transigen-shape-round1.md) ·
[Width extension](../../../docs/adr/0095-transigen-shape-encoder-width-extension.md) ·
[Registered shapes](../../util/training/shape_roster.py)

All fits use the same three frozen tenv cold-drug folds, the ordered
Denoise-DEG 5k axis, seed 42, frozen KPGT, numeric dose standardized on the
training fold, and measured same-plate DMSO. The Step 05.2 objective uses unit
expression and Δ MSE weights, Δ L1 weight 1.75, and weight 0.1 on the sum of
three native KL terms. The earliest minimum validation objective selects the
checkpoint, with patience 30 and a 1,000-epoch cap.
[Training release](../../03_data_train/03.7_tenv_denoise-5k/README.md) ·
[Step 05.2 protocol](../../05_tsg_loss/05.2_mse_l1/README.md)

## Formal results

The first matrix completed 24 fits with exact independent replay. Its three
`original` fits matched the Step 05.2 weight-1.75 learning curves, selected
checkpoints, validation and held-out predictions, Δ arrays, and all 27 fold
scores and eligibility counts. They selected epochs 129, 116, and 85 for
folds 0, 1, and 2. The extension reused those 24 sealed fits and independently
replayed six new fits at widths 2,000 and 4,000. The bounded report check passed
between the initial extension candidate and the published matrix. The repeat
check matched all 270 scores as well as the six new fits' selected epochs and
learning curves. [Initial receipt](log/347-formal.json) ·
[Combined receipt](log/350-formal.json) · [Repeat check](log/350-formal-repeat.json)

The 4,000-unit encoder had the highest mean Δ PCC (0.433229 vs. 0.431920 for
`original`) and lowest mean Δ RMSE (0.248576 vs. 0.249215). The 2,000- and
6,000-unit encoders and all three wider latents had lower mean Δ PCC and
higher mean Δ RMSE than `original`. All three wider encoders had higher mean
Δ E-distance and Δ MMD. Adding a second layer to both encoder and decoder gave
the lowest mean Δ E-distance (0.020687 vs. 0.025826) and Δ MMD (0.134288
vs. 0.154295). These are equal-fold means from one training seed.
[270 fold scores and eligibility counts](figure/data/fold-scores.csv) ·
[90 aggregates and paired differences](figure/data/delta-metrics.csv) ·
[Epoch objectives](figure/data/epoch-objectives.csv) · [Run ledger](log/runs.csv) ·
[KL diagnostics](log/kl-diagnostics.csv)

Step 06.2 reused the three `original` fits as the paired reference and kept
those expression dimensions for its prediction-architecture comparison.
[Step 06.2 results](../06.2_round2/README.md#formal-results) ·
[Step 06.2 design](../../../docs/adr/0096-transigen-shape-round2.md)

## Verify the accepted #347 source

The sealed matrix, audit, and report are under the installed
`output/formal/prior-347/`. The historical run names `widen_z-100`, `widen_enc`,
and `deepen-both` correspond to `original`, `widen_enc-6k`, and `deepen_both`
in the combined report. The sealed receipt records the training revision for
exact CUDA replay.
[Verification and replay guide](AGENTS.md#verify-the-accepted-347-source) ·
[Sealed receipt](log/347-formal.json) ·
[#347 release](https://github.com/Moke110/carps/tree/9dee5fd90ea6d2bed775655b0205875c8812938b/stage2/06_tsg_shape/06.1_round1)

## Reproduce the bounded smoke

The smoke covers `original` and the two added widths for two epochs on fold 0
with 64 train, 24 validation, and 24 test conditions.
The bounded report independently verified all 27 scores for the three fits
and reproduced both score figures byte for byte.
[Bounded smoke guide](AGENTS.md#reproduce-the-bounded-smoke) ·
[CUDA smoke receipt](log/349-smoke.json) · [Bounded report receipt](log/349-smoke-report.json) ·
[Pinned smoke recipe](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.1_round1/README.md#reproduce-the-bounded-smoke)

## Reproduce the formal matrix

The installed extension matrix, phase audits, and combined scored report
are under `output/formal/{matrix,matrix-audit,report}`. The pinned release
recipe retains the publisher and renderer used for the scientific receipt.
[Formal matrix guide](AGENTS.md#reproduce-the-formal-matrix) ·
[Combined receipt](log/350-formal.json) ·
[Pinned release recipe](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.1_round1/README.md#reproduce-the-formal-matrix)

## Replay the published figures

The current figures show the train and validation objectives for each fit,
held-out fold scores, and signed differences from `original` on the same fold
and seed. Bars show means, error bars show sample SD, and dots show individual
folds. The score panels include wJAC-2000 alongside JAC-top200 and the six Δ
metrics. The figure receipt binds the weighted score to the accepted
held-out predictions and matched controls.
[Learning curves](figure/learning-curve.png) · [Held-out scores](figure/metrics.png) ·
[Paired differences](figure/metrics-paired.png) · [Comparison protocol](../../../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md) ·
[Figure receipt](log/figure-render.json)

The renderer reads the three plotted CSVs. Its environment and the accepted
file digests are in the figure receipt.
[Figure replay guide](AGENTS.md#replay-the-published-figures) ·
[Epoch traces](figure/data/learning-curve.csv) · [Held-out values](figure/data/metrics.csv) ·
[Paired values](figure/data/metrics-paired.csv)
