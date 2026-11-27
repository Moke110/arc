# Step 10 CarPS-Transcriptome dose inputs

Step 10 tests whether supplying dose improves prediction for drug parents and
doses held out together. The successor crosses three historical held-parent
groups with 0.2, 1 and 10 µM on the ordered Denoise-DEG 5k axis. Each fold tests
the intersection, retains the historical development assignment at the other
two doses, and excludes conditions that meet exactly one holdout. Every
eligible condition is tested once across all nine folds.
[Shared tenv evaluation](../util/README.md#shared-tenv-evaluation) · [Rerun design](../../docs/adr/0108-cpst-loss-and-shared-evaluation.md)

| Fold | Held-parent group | Held dose (µM) | Train | Validation | Test | Excluded |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0 | 0.2 | 552 | 62 | 154 | 615 |
| 1 | 0 | 1 | 552 | 62 | 154 | 615 |
| 2 | 0 | 10 | 552 | 62 | 154 | 615 |
| 3 | 1 | 0.2 | 554 | 62 | 153 | 614 |
| 4 | 1 | 1 | 554 | 62 | 153 | 614 |
| 5 | 1 | 10 | 554 | 62 | 153 | 614 |
| 6 | 2 | 0.2 | 552 | 62 | 154 | 615 |
| 7 | 2 | 1 | 552 | 62 | 154 | 615 |
| 8 | 2 | 10 | 552 | 62 | 154 | 615 |

All nine CarPS-Transcriptome arms use frozen UniMol, the `XΔ-XΔ` backbone,
and the same objective. The eight dose-input arms apply scaling, FiLM or
concatenation on either side of the chemical projection. The TranSiGen
reference uses frozen KPGT with the author's `XX-XX` architecture and
objective. Same-plate measured DMSO defines both measured and predicted Δ.
Each fold's MMD bandwidth uses its complete training partition.
[UniMol decision](../../docs/adr/0100-carps-vae-unimol.md) · [MMD calibration](log/bandwidths.csv)

## Bounded CUDA verification

All ten arms passed the two-epoch CUDA gate on crossed fold 1, which holds
out parent group 0 at 1 µM, with seed 17. The formal matrix used seeds 17, 42
and 2026 and selected each checkpoint by minimum validation Loss.
[Bounded acceptance](log/405-bounded.json)

## Formal results

The formal matrix exactly replayed 270 checkpoints in independent
processes. The independent audit checked 2,430 score rows, 41,490
held-out condition rows and all 9 bandwidths. Values are means ± sample SD over
27 fold-and-seed coordinates per arm. The paired columns subtract the matched
`no-dose` result. The TranSiGen reference has absolute scores.
[Formal acceptance](log/405-formal.json) · [Fit lineage](log/runs.csv) · [Fold and seed scores](log/fold-scores.csv) · [Foldwise summaries](log/foldwise-summary.csv) · [All metrics](log/summary.csv)

The initial matrix scope stopped after 225 exact replays. The remaining
coordinates resumed on the same registered inputs and pinned runtime.
Systemd records bracket a 128.537-second interval between scopes.
Independent observations verify ten-second sampling in both active windows.
The initial scope reached 1.695 GiB in the sampled observations
before it stopped. The completed scope recorded a kernel memory peak of
1.737 GiB. The installed matrix audit retains the original
and resumed sampling and systemd records.
[Resource audit](log/405-formal.json)

| Arm | Δ Spearman ↑ | JAC-top200 ↑ | Δ PCC ↑ | Δ RMSE ↓ | Paired Δ PCC | Paired Δ RMSE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `TranSiGen` | 0.3683 ± 0.0181 | 0.1286 ± 0.0066 | 0.3722 ± 0.0222 | 0.2725 ± 0.0097 | — | — |
| `no-dose` | 0.3677 ± 0.0542 | 0.1463 ± 0.0214 | **0.3886 ± 0.0537** | **0.2454 ± 0.0076** | +0.0000 ± 0.0000 | +0.0000 ± 0.0000 |
| `Xmol-scaling` | 0.3115 ± 0.1150 | 0.1351 ± 0.0399 | 0.3320 ± 0.1197 | 0.2817 ± 0.0712 | -0.0566 ± 0.0989 | +0.0363 ± 0.0703 |
| `Xmol-FiLM` | 0.3257 ± 0.0630 | 0.1347 ± 0.0255 | 0.3433 ± 0.0635 | 0.2539 ± 0.0104 | -0.0453 ± 0.0473 | +0.0084 ± 0.0073 |
| `Zmol-scaling` | 0.3297 ± 0.1375 | 0.1430 ± 0.0443 | 0.3504 ± 0.1456 | 0.2626 ± 0.0442 | -0.0382 ± 0.1243 | +0.0172 ± 0.0438 |
| `Zmol-FiLM` | 0.3230 ± 0.0454 | 0.1351 ± 0.0201 | 0.3435 ± 0.0456 | 0.2530 ± 0.0072 | -0.0450 ± 0.0316 | +0.0076 ± 0.0062 |
| `Xdose-num` | 0.3204 ± 0.0712 | 0.1328 ± 0.0283 | 0.3380 ± 0.0730 | 0.2551 ± 0.0119 | -0.0506 ± 0.0682 | +0.0097 ± 0.0113 |
| `Xdose-MLP` | 0.3228 ± 0.0630 | 0.1346 ± 0.0261 | 0.3407 ± 0.0634 | 0.2535 ± 0.0096 | -0.0479 ± 0.0541 | +0.0081 ± 0.0080 |
| `Zdose-num` | 0.3153 ± 0.0693 | 0.1276 ± 0.0264 | 0.3343 ± 0.0713 | 0.2567 ± 0.0154 | -0.0542 ± 0.0637 | +0.0112 ± 0.0161 |
| `Zdose-MLP` | 0.3357 ± 0.0599 | 0.1347 ± 0.0238 | 0.3550 ± 0.0611 | 0.2525 ± 0.0122 | -0.0335 ± 0.0485 | +0.0071 ± 0.0123 |

`no-dose` had the highest mean Δ PCC and the lowest mean Δ RMSE across all
arms. `Zdose-MLP` led the dose-input arms on both metrics. The `no-dose`
configuration remains selected for subsequent CarPS-Transcriptome work.
Dose defines measured treatment conditions and the held-out partitions.
[Selection in ADR 0104](../../docs/adr/0104-carps-vae-no-dose.md)

## Figures

Each figure passed two independent CSV-only renders. Learning curves show
fitted Loss, and metric panels show fold-and-seed replicates with their means
and sample SD.
[Figure receipt](log/figure-render.json)

| Figure | Plot data |
| --- | --- |
| [Learning curves](figure/learning-curve.png) | [Training and validation Loss by epoch](figure/data/learning-curve.csv) |
| [Absolute metrics](figure/metrics.png) | [Scores by arm, fold and seed](figure/data/metrics.csv) |
| [Paired metrics](figure/metrics-paired.png) | [Differences from `no-dose`](figure/data/metrics-paired.csv) |

## Reproduce and verify

The formal runtime source is pinned at `20533e261485622436e5675bbb45250475d28eb6`. Use the
`output/formal-crossed-v1` bundle or pass its explicit materialized location
via `--bundle`. The predecessor published its three-fold study at a separate
pinned revision.
[Successor verification](AGENTS.md#verify-the-installed-release) · [Release checksum ledger](../util/training/cpst-receipts.sha256) · [Historical study](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/10_cpst_dose/README.md)
