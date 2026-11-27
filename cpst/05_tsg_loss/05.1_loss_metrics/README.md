# Step 05.1 — TranSiGen Delta loss metrics

Step 05.1 tested whether changing the form or weight of TranSiGen's Delta loss
improved held-out prediction. The 12 configurations crossed MSE, cosine, Pearson
correlation (PCC) and L1 with weights of 0.75, 1 and 1.25 on the
[tenv Denoise-DEG 5k release](../../03_data_train/03.7_tenv_denoise-5k/README.md).
Each configuration used the three frozen cold-drug folds, seed 42, frozen KPGT
vectors and the [Step 04 matched-plate protocol](../../04_tsg_axis/README.md).
Predicted and observed Delta both subtract the measured same-plate DMSO control,
which the model also receives as input.
[ADR 0089](../../../docs/adr/0089-transigen-delta-metric-weight.md) ·
[Spec #328](https://github.com/Moke110/carps/issues/328)

The objective was `0.1 * K + MSE_X1 + MSE_X2 + alpha * metric_Delta`, where
`K` sums the three native KL terms. MSE and L1 sum errors across conditions and
genes. The cosine and PCC losses sum per-condition dissimilarities and multiply
by the axis width. PCC centers each Delta vector over genes before normalization.
Both similarity losses clamp the L2 norms to a minimum of `1e-8`.
Each fit selected the checkpoint with minimum validation loss on its own objective,
with 30-epoch patience and a 1,000-epoch cap.
[Shared objective implementation](../../util/training/delta_objective.py)

## Accepted results and follow-up

| Delta loss | Weight | Mean Δ PCC ↑ | Mean Δ RMSE ↓ | Mean Δ E-distance ↓ | Mean Δ MMD ↓ |
| --- | ---: | ---: | ---: | ---: | ---: |
| Native MSE | 1 | 0.406610 | 0.259575 | 0.031911 | 0.179263 |
| L1 | 1.25 | 0.423530 | 0.254063 | 0.029019 | 0.166244 |

L1 at weight 1.25 led this grid on all four metrics shown. The means weight
the three held-out folds equally. [Complete equal-fold scores](figure/data/delta-metrics.csv) ·
[Fold scores and eligibility](log/fold-scores.csv) · [Frozen bandwidths](log/bandwidths.csv)

The accepted matrix contains 36 independently replayed fits and 324 fold scores.
All three MSE fits at weight 1 reproduced Step 04 exactly on selected and stopping
epochs, native objective traces, validation and test predictions, signed Delta
arrays, and all 27 reference fold scores and eligibility counts. Those fits served
as the reference for [Step 05.2](../05.2_mse_l1/README.md), which kept native Delta
MSE and added Delta L1. The [Step 05 decision](../README.md#selected-objective)
records the selected additive objective.
[Formal receipt and manifest hashes](log/330-formal.json) · [Run ledger](log/runs.csv)

The current [score figure](figure/metrics.png) and [paired differences](figure/metrics-paired.png)
show eight metrics. Differences subtract the native MSE score at weight 1 within the
same fold and seed. The figure update computed wJAC-2000 from the sealed predictions and
matched controls and verified parity on the sealed JAC-top200 scores. The
[learning curves](figure/learning-curve.png) show training and validation for every fit
through its stopping epoch, with 4,920 epochs in the [objective table](figure/data/epoch-objectives.csv).
[Render receipt](log/figure-render.json) · [CSV-only verification](AGENTS.md#verify-the-current-figures)

## Run a bounded comparison

The smoke study covered four two-epoch fits on fold 0 at weight 1.
[Smoke receipt](log/329-smoke.json) · [Smoke procedure](AGENTS.md#run-a-bounded-comparison)

## Reproduce the original formal matrix

The formal receipt pins the materialized `output/formal/matrix` and
`output/formal/report` manifests.
[Formal receipt](log/330-formal.json) ·
[Formal matrix procedure](AGENTS.md#reproduce-the-original-formal-matrix) ·
[Sealed replay procedure](AGENTS.md#replay-a-sealed-fit) ·
[Original publication commands](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/05_tsg_loss/05.1_loss_metrics/README.md#run-the-formal-matrix)
