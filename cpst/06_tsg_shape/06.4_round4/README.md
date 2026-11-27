# Step 06.4 TranSiGen measured encoder targets

Round 4 tested whether encoding measured signed Δ changes prediction
quality and whether the signed Δ decoder can be shared with reconstruction.
It retained the [Round 3 cold-drug protocol](../06.3_round3/README.md) on the
[tenv Denoise-DEG 5k release](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json).
All arms used measured same-plate DMSO, frozen KPGT and log dose standardized
on the training fold. The original prediction head uses the control latent
and drug–dose features to parameterize the Gaussian that supplies the
predicted latent. [Round 4 design and objectives](../../../docs/adr/0098-transigen-shape-round4.md)

| Arm | Measured encoders | Prediction readout |
| --- | --- | --- |
| `XX-XX` | Control and treated expression | Treated-expression decoder |
| `XX-XXΔ` | Control and treated expression | Independent signed Δ decoder |
| `X-XΔ` | Control only | Signed Δ decoder |
| `XΔ-XΔ` | Control and measured signed Δ | Signed Δ decoder shared with reconstruction |
| `XΔ-XΔΔ` | Control and measured signed Δ | Independent signed Δ decoder |

Direct Δ arms add decoded Δ to measured control to predict treated
expression. Both measured Δ arms reconstruct `X2 − X1` and align the
predicted Gaussian with its measured Δ Gaussian. The control-only
arm uses the control reconstruction and Gaussian prior KL, double
predicted Δ MSE and 1.75 times predicted Δ L1. The trace fields
`treatment_mse` and `treatment_kl` hold the corresponding measured Δ terms
in the two measured Δ arms and zeros in the control-only arm.

## Bounded CUDA gate

The two-epoch gate trained and independently replayed all 20 fits across
the five arms, folds 0–1 and training seeds 17 and 42. Each fold used
64 training, 24 validation and 24 held-out conditions. The independent
rescore matched all 180 metric values. The gate also checked the objective,
source lineage, CSV-only figure renders and resource use. [Gate receipt](log/358-gate.json)

The [Step 06.4 guide](AGENTS.md#bounded-cuda-gate) gives the bounded CUDA
procedure. The historical gate and renderer hashes use the
[original fit revision](https://github.com/Moke110/carps/blob/f387f2c370eb87cd27c02502aa98f165851491a7/stage2/06_tsg_shape/06.4_round4/README.md#bounded-cuda-gate).

## Formal results

The accepted comparison covers 45 coordinates across cold-drug folds
0–2 and training seeds 17, 42 and 2026. All 27 new fits passed separate-process
exact replay of the validation objective, validation predictions and
held-out predictions. The other 18 fits were verified against the
Step 06.3 release and its Step 06.1–06.2 lineage. Independent rescoring
matched all 405 held-out metric values. The fixed validation and test
sampling seeds were 22617 and 100042. [Formal acceptance receipt](log/359-formal.json)

| Arm | Δ PCC mean ± SD | Δ RMSE mean ± SD |
| --- | ---: | ---: |
| `XX-XX` | 0.426340 ± 0.014345 | 0.250777 ± 0.002121 |
| `XX-XXΔ` | 0.427947 ± 0.019165 | 0.241382 ± 0.002938 |
| `X-XΔ` | 0.347833 ± 0.102599 | 0.248650 ± 0.007588 |
| `XΔ-XΔ` | 0.428830 ± 0.018574 | **0.240752 ± 0.002598** |
| `XΔ-XΔΔ` | **0.429835 ± 0.017268** | 0.240819 ± 0.003052 |

Each arm mean weights its nine fold and seed scores equally. SD is the
sample standard deviation. The prespecified contrasts subtract scores
within each fold and seed before averaging. Each of the three held-out
folds has three training repeats.

| Paired contrast | Δ PCC | Δ RMSE | Δ Spearman | Δ E-distance |
| --- | ---: | ---: | ---: | ---: |
| `XΔ-XΔΔ − XX-XXΔ` | +0.001888 ± 0.007564 | −0.000563 ± 0.001480 | +0.001995 ± 0.007217 | +0.001090 ± 0.001945 |
| `XΔ-XΔ − XΔ-XΔΔ` | −0.001005 ± 0.010334 | −0.000067 ± 0.000836 | −0.000546 ± 0.010250 | −0.000210 ± 0.002162 |

For both contrasts, the paired sample SD exceeds the absolute mean
difference on all four metrics. The independent measured Δ decoder had
the highest mean Δ PCC. The shared decoder had the lowest mean Δ RMSE
and became the fixed CarPS-Transcriptome alpha under ADR 0099. Step 07
verifies that the release reproduces the nine shared-decoder fits.
[Architecture decision](../../../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md) ·
[Alpha release and checkpoint parity](../../07_cpst_alpha/README.md#checkpoint-and-replay-evidence)

- [Nine metric summaries](figure/data/delta-metrics.csv) and
  [405 held-out scores with eligibility counts](figure/data/fold-scores.csv)
- [72 paired values](figure/data/primary-contrasts.csv),
  [24 foldwise summaries](figure/data/primary-contrast-by-fold.csv) and
  [eight overall summaries](figure/data/primary-contrast-summary.csv)
- [Fit ledger](log/runs.csv), [epoch objectives](figure/data/epoch-objectives.csv)
  and [fold bandwidths](log/bandwidths.csv)
- [Learning curves](figure/learning-curve.png), [held-out metrics](figure/metrics.png),
  [differences from `XX-XX`](figure/metrics-paired.png) and
  [original paired contrast figure](https://github.com/Moke110/carps/blob/09da1ead4ee31d5c8d84ac60d252a4a2f3cc209f/stage2/06_tsg_shape/06.4_round4/figure/primary-contrasts.png)
- [Formal matrix resource log](log/359-resources.jsonl) and
  [coordinator resource log](log/359-coordinator-resources.jsonl)

## Reproduce the formal study

The [formal receipt](log/359-formal.json) pins the matrix and report manifests
under `output/formal/` in a materialized checkout. The
[Step 06.4 guide](AGENTS.md#reproduce-the-formal-study) covers fitting and
verification, including reference lineage and independent rescoring. The
[original formal publication and verification commands](https://github.com/Moke110/carps/blob/09da1ead4ee31d5c8d84ac60d252a4a2f3cc209f/stage2/06_tsg_shape/06.4_round4/README.md#reproduce-the-formal-study)
preserve the original layout, paired table audits and figure renders.

## Replay a sealed Step 06.4 fit

The formal receipt pins the matrix at `f387f2c370eb87cd27c02502aa98f165851491a7`.
The registrations bind the original `stage2/` paths and code digests. The
[Step 06.4 guide](AGENTS.md#replay-a-sealed-step-064-fit) shows how to
replay the sealed shared-decoder fit for fold 0, seed 17 from a scratch copy.

## Replay the published figures

The current figures draw from the [epoch traces](figure/data/learning-curve.csv),
[held-out scores](figure/data/metrics.csv) and
[paired differences](figure/data/metrics-paired.csv). Paired plots subtract
`XX-XX` within each fold and seed. The [figure receipt](log/figure-render.json)
pins the accepted fit evidence, plotted CSVs, PNGs and renderer versions.
The [Step 06.4 guide](AGENTS.md#replay-the-published-figures) gives the
CSV-only replay for these three figures.
