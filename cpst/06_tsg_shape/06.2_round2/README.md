# Step 06.2 TranSiGen prediction architectures

Step 06.2 tests eight ways to predict treated expression from a measured
control and a chemical perturbation. It keeps the Step 06.1 original
expression architecture (1,200 encoder units, 800 decoder units, and
100-dimensional latents) and compares changes to the readout, prediction
head, or encoder sharing. All arms use three frozen tenv cold-drug folds,
the ordered Denoise-DEG 5k axis, seed 42, measured same-plate DMSO, frozen
KPGT, and 32 repeats of the log-dose standardized on each training fold.
[Study design](../../../docs/adr/0096-transigen-shape-round2.md) ·
[Step 06.1 baseline](../06.1_round1/README.md#formal-results) ·
[Training release](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json)

| Architecture | Registered change |
| --- | --- |
| `original` | Step 06.1 original architecture and Step 05.2 objective. |
| `output-delta` | Add an independent decoder with a signed linear Δ output from the predicted treated latent. |
| `output-delta-nox2` | Repurpose the treated decoder for a signed Δ output and drop the treated reconstruction MSE. |
| `z-interact_MLP` | Predict Gaussian mean and log variance with a 532→145→200 ReLU MLP. |
| `z-interact_FiLM` | Use a 432→136→200 ReLU conditioner with identity-initialized modulation and a 100→200 Gaussian projection. |
| `shared-z` | Share encoder trunks and their BatchNorm state, with separate Gaussian projections. |
| `z-res` | Add the sampled control latent to the predicted Gaussian mean. |
| `shared-z-res` | Combine the shared encoder trunk and the prediction-mean residual. |

Both signed Δ arms predict treated expression by adding the signed output
to the measured control. The `output-delta-nox2` objective keeps the measured
treated encoder and all three native KL terms and uses
`0.1 K + MSE_X1 + 2 MSE_Δ + 1.75 L1_Δ` for fitting and checkpoint selection.
The other arms keep the Step 05.2 objective with unit weights on control,
treated, and Δ MSE. Round 2 FiLM has no normalization after modulation.
[Step 05.2 objective](../../05_tsg_loss/05.2_mse_l1/README.md) ·
[Architecture implementation](../../util/training/architecture_model.py)

## Formal results

FiLM had the highest mean Δ PCC and improved both primary metrics on two
folds. The independent signed Δ decoder had the lowest mean Δ RMSE.
Both signed Δ arms reduced Δ RMSE on every fold and had lower mean Δ PCC
than `original`. The prediction-mean residual also reduced Δ RMSE on all
three folds and raised mean Δ PCC. The residual changes the parameterization
of the existing linear mean head.

The means give equal weight to the three folds at seed 42. Each paired
difference subtracts `original` on the same fold and is shown as mean ±
sample SD. Higher Δ PCC and lower Δ RMSE are favorable.

| Architecture | Mean Δ PCC | Mean Δ RMSE | Paired Δ PCC difference | Paired Δ RMSE difference |
| --- | ---: | ---: | ---: | ---: |
| `original` | 0.431920 | 0.249215 | +0.000000 ± 0.000000 | +0.000000 ± 0.000000 |
| `output-delta` | 0.423164 | 0.241952 | −0.008756 ± 0.025504 | −0.007263 ± 0.005194 |
| `output-delta-nox2` | 0.415948 | 0.241974 | −0.015972 ± 0.010496 | −0.007241 ± 0.003617 |
| `z-interact_MLP` | 0.423109 | 0.252523 | −0.008812 ± 0.007767 | +0.003308 ± 0.003277 |
| `z-interact_FiLM` | 0.439568 | 0.247354 | +0.007647 ± 0.013368 | −0.001861 ± 0.005067 |
| `shared-z` | 0.408751 | 0.254451 | −0.023169 ± 0.014882 | +0.005236 ± 0.002960 |
| `z-res` | 0.436771 | 0.246442 | +0.004851 ± 0.013173 | −0.002773 ± 0.003485 |
| `shared-z-res` | 0.424692 | 0.247712 | −0.007228 ± 0.009384 | −0.001503 ± 0.003217 |

All 21 new fits passed independent exact replay of validation predictions,
validation objective, and held-out predictions. The three sealed Step 06.1
originals reproduced all 27 reference scores and eligibility counts.
Independent scoring reproduced all 216 fold scores and 72 aggregates over
the nine formal metrics.
[Formal receipt](log/353-formal.json) ·
[Fold scores and eligibility](figure/data/fold-scores.csv) ·
[Aggregates and paired differences](figure/data/delta-metrics.csv) ·
[Epoch objectives](figure/data/epoch-objectives.csv) ·
[Fit ledger](log/runs.csv)

Step 06.3 carried both signed readouts and FiLM into a four-seed
combination study. Its accepted run reused six seed-42 signed-readout fits
from Round 2 and trained the FiLM arms with parameter-free LayerNorm.
Step 06.4 took the signed Δ readout forward with the original Gaussian
prediction head and supplied the architecture chosen for CarPS-Transcriptome.
[Combination study](../06.3_round3/README.md#formal-results) ·
[CarPS-Transcriptome selection](../../README.md#selected-carps-transcriptome-configuration)

## Reproduce the bounded CUDA smoke

The accepted smoke used all eight arms on fold 0 for two epochs, with 64
training, 24 validation, and 24 test conditions. All eight runs exactly
replayed the validation objective and validation and test predictions.
The smoke also independently verified 72 scores and rendered three
figures twice from their CSVs.
[Smoke receipt](log/352-smoke.json) ·
[Smoke commands](AGENTS.md#reproduce-the-bounded-cuda-smoke)

## Reproduce the formal matrix

The installed release keeps the matrix, phase audits, and scored report
under `output/formal/{matrix,matrix-audit,report}`. The original #353
publication preserves the `stage2/` source layout, the formal verifier, and
the three figures it expects. Formal verification uses the installed products,
the frozen training and calibration inputs, and the verified Step 06.1
predecessor matrix. The #353 receipt records the fit revision and resource checks.
[Formal verification guide](AGENTS.md#reproduce-the-formal-matrix) ·
[Pinned formal recipe](https://github.com/Moke110/carps/blob/87564cda35c7887de2cf6c79d1b2f0221fc368e8/stage2/06_tsg_shape/06.2_round2/README.md#reproduce-the-formal-matrix) ·
[Formal receipt](log/353-formal.json) ·
[Step 06.1 predecessor](../06.1_round1/README.md#verify-the-accepted-347-source) ·
[Installed release inventory](../../util/training/verification.md)

## Replay the published figures

The current figures show training and validation objectives, held-out
scores, and paired differences from `original` within fold and seed.
The score and difference panels use fold dots, means, and sample SD.
They include wJAC-2000 from the accepted predictions and matched controls
alongside JAC-top200 and the six Δ metrics.
[Learning curves](figure/learning-curve.png) ·
[Held-out scores](figure/metrics.png) ·
[Paired differences](figure/metrics-paired.png) ·
[Figure receipt](log/figure-render.json) ·
[CSV-only replay commands](AGENTS.md#replay-the-published-figures)

[Epoch traces](figure/data/learning-curve.csv) ·
[Held-out values](figure/data/metrics.csv) ·
[Paired values](figure/data/metrics-paired.csv)
