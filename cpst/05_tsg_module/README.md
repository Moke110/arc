# Step 05 TranSiGen modules

This study compares five modules on the same tenv Denoise-DEG 5k training
release, measured same-plate DMSO, frozen KPGT and fold-local standardized
log dose. Each module uses the five shared cold-drug folds and training
seeds 17, 42, 43, 2026 and 2027, for 125 new coordinates.
[ADR 0110](../../docs/adr/0110-cpst-transigen-module-rerun.md) ·
[Training and evaluation commands](AGENTS.md)

| Module | Measured encoders | Prediction |
| --- | --- | --- |
| `XX-XX` | Control and treated expression | Treated expression |
| `XX-XXΔ` | Control and treated expression | Independent signed Δ decoder |
| `X-XΔ` | Control only | Signed Δ decoder |
| `XΔ-XΔ` | Control and measured signed Δ | Shared signed Δ decoder |
| `XΔ-XΔΔ` | Control and measured signed Δ | Independent signed Δ decoder |

Training and checkpoint selection minimize the sum of control reconstruction
MSE, signed Δ reconstruction MSE, predicted Δ MSE and 0.1 times the KL sum.
The control-only module uses the control prior KL, control reconstruction and
predicted Δ MSE. The other modules use both prior KL terms and a
predicted-to-measured Gaussian alignment KL. For the two treated-expression
reconstruction modules, subtract the same measured control from reconstructed
and observed treated expression. Their squared reconstruction error therefore
equals the treated-expression squared error algebraically.

The comparison figure subtracts each metric's mean over all valid displayed
`XX-XX` coordinates from every displayed fold-and-seed score. The reference
group has its own centered dots and sample SD. The independent primary
contrasts still subtract scores within the same fold and seed.

## Bounded validation

All five modules passed their two-epoch CUDA fit and exact fresh-process replay
on the complete first fold with seed 17. An independent calculation reproduced
the selected checkpoints' full validation objectives. The bounded report
independently reproduced 35 scores and the paired tables. Two CSV-only renders
reproduced all three figures. The formal register-only run verified 125 new
registered coordinates. [Validation receipt](log/410-gate.json)

## Formal results

All 125 fresh fits passed independent CUDA replay of the selected checkpoint's
validation objective, validation predictions and held-out predictions. A
separate scoring process reproduced all 875 held-out scores from the sealed
predictions. The report retains 34,575 condition-level rows. Each module
summary gives the mean and sample SD of its 25 fold-and-seed scores.

| Module | Δ Spearman | JAC-top200 | Δ PCC | Δ RMSE |
| --- | ---: | ---: | ---: | ---: |
| `XX-XX` | **0.39111 ± 0.02161** | 0.14113 ± 0.00741 | **0.40125 ± 0.01631** | 0.25802 ± 0.00476 |
| `XX-XXΔ` | 0.36452 ± 0.02780 | 0.14483 ± 0.00775 | 0.38817 ± 0.02313 | 0.24801 ± 0.00638 |
| `X-XΔ` | 0.36265 ± 0.02681 | **0.14666 ± 0.00758** | 0.38615 ± 0.02257 | **0.24723 ± 0.00529** |
| `XΔ-XΔ` | 0.36929 ± 0.02306 | 0.14588 ± 0.00671 | 0.39182 ± 0.01802 | 0.24763 ± 0.00646 |
| `XΔ-XΔΔ` | 0.36713 ± 0.02438 | 0.14589 ± 0.00665 | 0.39071 ± 0.01959 | 0.24744 ± 0.00615 |

`XX-XX` had the highest mean Δ Spearman and Δ PCC. The control-only
`X-XΔ` had the highest mean signed Jaccard and the lowest mean Δ RMSE.
Both primary matched contrasts had sample SDs larger than the absolute mean
difference on Δ Spearman, Δ PCC, and Δ RMSE. The shared measured-Δ decoder
had higher Δ E-distance than the independent measured-Δ decoder in the
same-coordinate contrast, by 0.002879 on average with a sample SD of 0.001762.

- [Seven metric summaries](figure/data/delta-metrics.csv) and
  [875 held-out scores with eligibility counts](figure/data/fold-scores.csv)
- [200 primary contrast scores](figure/data/primary-contrasts.csv),
  [40 foldwise summaries](figure/data/primary-contrast-by-fold.csv) and
  [eight overall summaries](figure/data/primary-contrast-summary.csv)
- [Learning curves](figure/learning-curve.png), [held-out metrics](figure/metrics.png)
  and [differences from the XX-XX mean](figure/metrics-paired.png)
- [Fit ledger](log/runs.csv), [fold calibration](log/bandwidths.csv),
  [formal acceptance](log/411-formal.json), [installed acceptance](log/411-installed.json)
  and [CSV replay receipt](log/figure-render.json)

Use the [installed replay procedure](AGENTS.md#accept-and-replay-the-installed-publication)
for independent scoring and the two CSV-only figure renders.

## Historical reports

The prior loss studies and the four module rounds remain available at the
revision before the Step 05 cutover. Use their pinned scientific runtime
and saved inputs for historical replays. The existing Step 07–10 model
decisions derive from those historical releases.
[Loss study](https://github.com/Moke110/carps/blob/325be5dd45e2b670f1e01882e225aa5debbeaab7/cpst/05_tsg_loss/README.md) ·
[Module rounds and width extension](https://github.com/Moke110/carps/blob/325be5dd45e2b670f1e01882e225aa5debbeaab7/cpst/06_tsg_module/README.md) ·
[Local evidence locations](../util/training/verification.md)
