# Step 10 CarPS-Transcriptome dose inputs

Step 10 tested whether supplying measured dose improved prediction for drug parents and doses held out together. It fixed the ordered Denoise-DEG 5k axis selected in Step 09 and used the frozen tenv `cold_dose_drug` folds. Within each fold, the test set contains conditions at the held dose for held parents. Training and validation contain other parents at the other two doses. Conditions meeting just one holdout are excluded. Each dose is tied to a different held-parent group, so differences across folds combine dose and parent effects. [Step 09 selection](../09_cpst_axis/README.md#formal-results) · [Frozen training release](../03_data_train/03.7_tenv_denoise-5k/README.md) · [Study design](../../docs/adr/0103-carps-vae-dose-inputs.md)

| Fold | Held dose (µM) | Train | Validation | Test | Excluded |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.2 | 552 | 62 | 154 | 615 |
| 1 | 1 | 554 | 62 | 153 | 614 |
| 2 | 10 | 552 | 62 | 154 | 615 |

All nine CarPS-Transcriptome arms used frozen UniMol, the `XΔ-XΔ` backbone with a shared signed Δ decoder, and the same objective. The eight dose-input arms applied scaling, FiLM or concatenation on either side of the chemical projection. The TranSiGen reference used frozen KPGT and the author's `XX-XX` architecture and objective without a dose input. Its absolute scores reflect the combined architecture, representation and objective differences. [Alpha backbone](../07_cpst_alpha/README.md) · [UniMol decision](../../docs/adr/0100-carps-vae-unimol.md)

The study selected checkpoints by minimum validation objective. It computed measured and predicted Δ by subtracting the measured same-plate DMSO profile from both treated profiles and used full-axis Δ PCC and Δ RMSE as primary metrics. Each fold's MMD bandwidth was calibrated on its complete training partition and shared across arms and seeds. [MMD calibration](log/bandwidths.csv)

## Bounded CUDA verification

All ten arms passed the two-epoch gate on fold 0 with seed 17 using the first 64 training, 24 validation and 24 test conditions. The receipt records the registered model interfaces, exact replay of the validation objective and validation and test predictions, report and figure verification, and resource measurements. Its configuration digest binds the 90-coordinate formal roster. [Bounded acceptance receipt](log/377-bounded.json)

## Formal results

**Step 10 selected `no-dose` for subsequent CarPS-Transcriptome work.** It had the highest mean full-axis Δ PCC and lowest mean Δ RMSE of the ten arms. `Zdose-MLP` led the dose-input arms on both metrics. The selected prediction head takes the control latent and projected frozen UniMol representation. Dose continues to identify measured treatment conditions and define the held-out partitions. [Selected configuration](../README.md#selected-carps-transcriptome-configuration) · [Decision in ADR 0104](../../docs/adr/0104-carps-vae-no-dose.md)

All 90 fits across three folds and seeds 17, 42 and 2026 exactly replayed their validation objective, validation predictions and held-out predictions in independent processes. A separate audit reconstructed the matched-control Δ scores and checked 900 score rows and 13,830 held-out condition rows, with absolute and relative numeric tolerances of 10⁻¹². Table entries give equal-weight mean ± sample SD over the nine fold and seed coordinates per arm. Paired differences subtract the `no-dose` score for the same fold and seed and apply to the CarPS-Transcriptome arms. [Formal acceptance receipt](log/378-formal.json) · [Fit lineage](log/runs.csv) · [Fold and seed scores](log/fold-scores.csv) · [Foldwise summaries](log/foldwise-summary.csv) · [All metrics](log/summary.csv)

| Arm | Δ PCC ↑ | Δ RMSE ↓ | Paired Δ PCC | Paired Δ RMSE |
| --- | ---: | ---: | ---: | ---: |
| TranSiGen | 0.3640 ± 0.0166 | 0.2767 ± 0.0037 | — | — |
| `no-dose` | **0.3913 ± 0.0235** | **0.2472 ± 0.0046** | +0.0000 ± 0.0000 | +0.0000 ± 0.0000 |
| `Xmol-scaling` | 0.3421 ± 0.1135 | 0.2736 ± 0.0479 | −0.0491 ± 0.0979 | +0.0264 ± 0.0503 |
| `Xmol-FiLM` | 0.3342 ± 0.0618 | 0.2569 ± 0.0068 | −0.0570 ± 0.0593 | +0.0097 ± 0.0050 |
| `Zmol-scaling` | 0.3392 ± 0.1476 | 0.2679 ± 0.0424 | −0.0521 ± 0.1303 | +0.0208 ± 0.0448 |
| `Zmol-FiLM` | 0.3418 ± 0.0414 | 0.2543 ± 0.0053 | −0.0494 ± 0.0395 | +0.0071 ± 0.0053 |
| `Xdose-num` | 0.3190 ± 0.0708 | 0.2606 ± 0.0080 | −0.0722 ± 0.0668 | +0.0135 ± 0.0086 |
| `Xdose-MLP` | 0.3313 ± 0.0635 | 0.2567 ± 0.0078 | −0.0599 ± 0.0653 | +0.0096 ± 0.0074 |
| `Zdose-num` | 0.3512 ± 0.0544 | 0.2561 ± 0.0081 | −0.0400 ± 0.0425 | +0.0089 ± 0.0108 |
| `Zdose-MLP` | 0.3700 ± 0.0324 | 0.2511 ± 0.0041 | −0.0212 ± 0.0254 | +0.0039 ± 0.0056 |

## Figures

Two independent CSV-only renders reproduced the three figures. The learning curves show each model's own objective, with transient spikes compressing some panels. The absolute and paired figures show fold and seed replicates, means and sample SD. [Figure receipt](log/figure-render.json)

| Figure | Plot data |
| --- | --- |
| [Learning curves](figure/learning-curve.png) | [Training and validation objectives by epoch](figure/data/learning-curve.csv) |
| [Absolute metrics](figure/metrics.png) | [Scores by arm, fold and seed](figure/data/metrics.csv) |
| [Paired metrics](figure/metrics-paired.png) | [Differences from `no-dose` by fold and seed](figure/data/metrics-paired.csv) |

## Reproduce and verify

The active figures and score tables retain the accepted values for the
current metrics. Use the pinned publication for full scientific replay.
[Sealed publication](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/10_cpst_dose/README.md)

The installed `output/formal` bundle holds the sealed fit matrix, audit evidence and full scored report. The formal receipt and shared checksum ledger bind its manifests to the original tables and figures. [Release checksum ledger](../util/training/cpst-receipts.sha256)

The operating guide covers fitting, report reconstruction and CSV-only figure replay. [Run the bounded CUDA gate](AGENTS.md#run-the-bounded-cuda-gate) · [Reproduce and verify the formal study](AGENTS.md#reproduce-and-verify-the-formal-study) · [Rebuild a report from sealed fits](AGENTS.md#rebuild-a-report-from-sealed-fits) · [Verify the installed release](AGENTS.md#verify-the-installed-release)
