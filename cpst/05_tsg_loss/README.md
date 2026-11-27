# TranSiGen loss studies

Step 05 tests which loss terms and weights improve held-out drug response
prediction. All five studies use the tenv Denoise-DEG 5k training release,
three frozen cold-drug folds, seed 42 and frozen KPGT vectors. Measured
same-plate DMSO is the model input `X1` and the common reference for
predicted and observed signed Δ. The results apply when that matched
control is available.
[Training release](../03_data_train/03.7_tenv_denoise-5k/README.md) ·
[Matched-control comparison](../04_tsg_axis/README.md)

## Loss studies

| Experiment | Accepted evidence and follow-up |
| --- | --- |
| [05.1 Delta loss forms](05.1_loss_metrics/README.md) | Crossed MSE, cosine, PCC and L1 with three weights. Native MSE at weight 1 exactly reproduced the Step 04 reference. L1 at 1.25 led the grid on mean Δ PCC (0.423530), Δ RMSE (0.254063), Δ E-distance and Δ MMD. Step 05.2 tested L1 alongside native MSE. |
| [05.2 Additive Delta L1](05.2_mse_l1/README.md) | Added Δ L1 to the native objective. Weight 1.75 raised mean Δ PCC from 0.406610 to 0.431920 and lowered mean Δ RMSE from 0.259575 to 0.249215 against the zero-weight reference. Weight 1.5 gave slightly lower Δ E-distance and Δ MMD. Step 05.3 fixed the L1 weight at 1.75. |
| [05.3 Reconstruction MSE](05.3_mse-l1_mse-x/README.md) | Varied one shared weight for both expression reconstruction terms. Unit weights exactly reproduced the Step 05.2 reference at L1 weight 1.75 and led on mean Δ PCC and Δ RMSE. Zero reconstruction weight had the lowest mean Δ E-distance (0.020303) and Δ MMD (0.143878). |
| [05.4 Bounded Delta PCC](05.4_mse-l1_pcc/README.md) | Added a PCC penalty to the unit-reconstruction objective. Zero added weight led on mean Δ PCC and Δ RMSE. |
| [05.5 Bounded Delta cosine](05.5_mse-l1_cos/README.md) | Added a cosine penalty to the same objective. Zero added weight led on mean Δ PCC and Δ RMSE. Weight 0.3 raised mean Jaccard overlap among the 100 most downregulated genes from 0.183801 to 0.185508. |

The formal receipts cover 126 independently replayed fits. Each similarity
sweep used fresh zero-weight fits that exactly reproduced the Step 05.3
unit-reconstruction reference.
[05.1 receipt](05.1_loss_metrics/log/330-formal.json) ·
[05.2 receipt](05.2_mse_l1/log/333-formal.json) ·
[05.3 receipt](05.3_mse-l1_mse-x/log/336-formal.json) ·
[05.4 receipt](05.4_mse-l1_pcc/log/339-formal.json) ·
[05.5 receipt](05.5_mse-l1_cos/log/339-formal.json)

## Selected objective

`L = 0.1 * K + MSE_X1 + MSE_X2 + MSE_Δ + 1.75 * L1_Δ`

`K` sums the three native KL terms. The selected Step 05.2 coordinate
keeps unit weights on control reconstruction, treated reconstruction and
Δ MSE. It attained the highest equal-fold mean Δ PCC and lowest equal-fold
mean Δ RMSE in Step 05, tied by the later reference fits.
[Step 05.2 equal-fold scores](05.2_mse_l1/figure/data/delta-metrics.csv)

Step 06.1 carried this objective into the architecture comparison and
exactly reproduced the three selected Step 05.2 fits with the original
architecture. The downstream CarPS-Transcriptome alpha kept the unit
predicted Δ MSE and 1.75 predicted Δ L1 weights with a shared signed Δ
decoder and a measured Δ reconstruction term.
[Step 06.1 architecture results](../06_tsg_shape/06.1_round1/README.md#formal-results) ·
[Fixed CarPS-Transcriptome alpha](../07_cpst_alpha/README.md#fixed-model-and-public-api)

## Verify the current figures

The five studies each publish learning curves, held-out scores and paired
differences. The comparison figures show seven metrics
and retain zero on every score axis. Bars show means, dots mark
individual fold and seed values, and error bars show sample SD. Paired
differences subtract each study's reference on the same fold and seed.
The receipts record two independent renders of all 15 PNGs from plotted CSVs.
[05.1 render](05.1_loss_metrics/log/figure-render.json) ·
[05.2 render](05.2_mse_l1/log/figure-render.json) ·
[05.3 render](05.3_mse-l1_mse-x/log/figure-render.json) ·
[05.4 render](05.4_mse-l1_pcc/log/figure-render.json) ·
[05.5 render](05.5_mse-l1_cos/log/figure-render.json) ·
[CSV-only verification guide](AGENTS.md#verify-the-current-figures)

## Rebuild the current comparisons from accepted fits

The builder checks accepted fit and report manifests, reconstructs signed Δ
from fixed predictions and matched controls, and checks JAC-top200 parity. It
uses complete training and validation traces to produce the comparison CSVs,
figures and receipts.
[Accepted-fit rebuild guide](AGENTS.md#rebuild-the-current-comparisons-from-accepted-fits) ·
[Original publication and figures](https://github.com/Moke110/carps/tree/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/05_tsg_loss)

## Reproduce the original formal Delta similarity sweeps

The PCC and cosine studies share a 36-fit formal matrix and have separate
scoring receipts. The Step 05 guide covers the joint procedure and each
study guide provides the pinned route for replaying an individual sealed fit.
[Joint formal sweep guide](AGENTS.md#reproduce-the-original-formal-delta-similarity-sweeps) ·
[Pinned PCC replay](05.4_mse-l1_pcc/AGENTS.md#replay-a-sealed-fit) ·
[Pinned cosine replay](05.5_mse-l1_cos/AGENTS.md#replay-a-sealed-fit)
