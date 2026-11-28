# Step 06.3 TranSiGen architecture combinations

Round 3 tested whether FiLM improves either signed Δ readout from Round 2.
The comparison used three frozen tenv cold-drug folds on the ordered
Denoise-DEG 5k axis, measured same-plate DMSO, frozen KPGT and 32 repeats of
log dose standardized on the training fold. Seeds 17, 42, 2026 and 43 repeated
training within each fold. [Round 3 design](../../../docs/adr/0097-transigen-shape-round3.md) ·
[Round 2 protocol](../../../docs/adr/0096-transigen-shape-round2.md)

| Architecture | Readout and prediction head |
| --- | --- |
| `0_original` | Original treated-expression readout and Gaussian prediction head |
| `1_output-delta` | Independent signed Δ decoder with the original prediction head |
| `2_output-delta-no-x2` | Treated decoder emits signed Δ with the original prediction head |
| `3_z-interact_FiLM` | Original treated-expression readout with normalized FiLM |
| `4_1+3` | Independent signed Δ decoder with normalized FiLM |
| `5_2+3` | Treated decoder emits signed Δ with normalized FiLM |

The four signed Δ arms predict treated expression as measured control plus
predicted Δ. Arms 2 and 5 keep the measured treated encoder and all three KL
terms and replace the treated reconstruction MSE with a second Δ MSE term.
The other arms use the Step 05.2 objective. [Loss selection](../../05_tsg_loss/05.2_mse_l1/README.md)

## First attempt and normalization

The initial bounded check passed all 24 independent fit replays. The first
formal attempt then stopped with prediction KL overflow in eight coordinates
from arms 4 and 5. It independently replayed 52 successful new fits before
ending without a formal publication. [Initial smoke](log/355-smoke.json) ·
[First attempt](log/356-blocked.json)

The revision kept identity initialization for FiLM modulation and added
`LayerNorm(100, elementwise_affine=False)` between modulation and the Gaussian
projection in arms 3–5. The revised two-epoch gate replayed 24 fits across all
six arms on folds 0 and 1 and seeds 17 and 42. Its objective audit checked 96
weighted objective rows. Formal execution paused after 13 sealed fits while
the revised gate ran, then resumed after checking those fits against their
pre-gate hashes. [Revised smoke](log/356-normalized-smoke.json) ·
[Objective audit](log/356-normalized-smoke-objective-audit.json) ·
[Formal resource timeline](log/356-coordinator-resources.jsonl)

## Formal results

All 36 normalized FiLM fits passed exact independent replay of the validation
objective, validation predictions and held-out predictions. This covers all
eight failed coordinates from the first attempt. The other 36 fits passed
source and artifact verification. They comprise nine seed-42 fits from Rounds
1 and 2 and 27 fits from the first Round 3 attempt. Independent rescoring
reproduced the 81 earlier seed-42 scores and their eligibility records and
covered all 648 scores in the revised comparison. [Formal receipt](log/356-formal.json) ·
[Fit ledger](log/runs.csv) · [Independent table audit](log/356-formal-table-audit.json)

Each arm has 12 equally weighted fold and training-seed coordinates. The paired
differences subtract `0_original` on the same fold and seed and are shown as
mean ± sample SD. The four seeds repeat fitting within the same held-out splits.

| Architecture | Mean Δ PCC | Mean Δ RMSE | Paired Δ PCC difference | Paired Δ RMSE difference |
| --- | ---: | ---: | ---: | ---: |
| `0_original` | 0.424154 | 0.251615 | +0.000000 ± 0.000000 | +0.000000 ± 0.000000 |
| `1_output-delta` | 0.429189 | 0.241190 | +0.005035 ± 0.016649 | −0.010425 ± 0.004631 |
| `2_output-delta-no-x2` | 0.421896 | 0.241106 | −0.002258 ± 0.018312 | −0.010509 ± 0.004267 |
| `3_z-interact_FiLM` | 0.396089 | 0.259210 | −0.028065 ± 0.017481 | +0.007595 ± 0.005446 |
| `4_1+3` | 0.428704 | 0.242163 | +0.004550 ± 0.013717 | −0.009452 ± 0.004250 |
| `5_2+3` | 0.400866 | 0.243244 | −0.023288 ± 0.016215 | −0.008371 ± 0.004231 |

The independent signed Δ decoder gave the highest mean Δ PCC and improved
mean Δ RMSE relative to the original. Adding normalized FiLM to that decoder
lowered mean Δ PCC by 0.000485 and raised mean Δ RMSE by 0.000973. Adding FiLM
to the other signed readout also lowered mean Δ PCC and raised mean Δ RMSE.
Both combinations improved the two means relative to FiLM alone. The component
and descriptive interaction contrasts use the same paired coordinates. [Contrast summaries](figure/data/primary-contrast-summary.csv) ·
[Coordinate contrasts](figure/data/primary-contrasts.csv)

Step 06.4 kept the original Gaussian prediction head and verified reuse of 18
fits from `0_original` and `1_output-delta` to compare measured encoder targets
and decoder sharing. [Round 4 design](../../../docs/adr/0098-transigen-shape-round4.md) ·
[Round 4 results](../06.4_round4/README.md#formal-results)

[Held-out scores](figure/metrics.png) · [Paired differences](figure/metrics-paired.png) ·
[Learning curves](figure/learning-curve.png) · [All nine sealed metric summaries](figure/data/delta-metrics.csv) ·
[Fold and seed scores](figure/data/fold-scores.csv)

## Reproduce the bounded CUDA smoke

The smoke uses the first 64 training, 24 validation and 24 held-out registered
conditions per fold. Validation and test sampling seeds are 22617 and 100042.
[Current smoke procedure](AGENTS.md#reproduce-the-bounded-cuda-smoke) ·
[Pinned original smoke procedure](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.3_round3/README.md#reproduce-the-bounded-cuda-smoke)

## Reproduce the revised formal run

The owning guide covers new registrations and verification of the installed
release with the pinned verifier. [Formal run and verification](AGENTS.md#reproduce-the-revised-formal-run) ·
[Original formal publication](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.3_round3/README.md#reproduce-the-revised-formal-run) ·
[Installed release inventory](../../util/training/verification.md) ·
[Shared TranSiGen APIs](../../util/README.md#transigen-training-and-evaluation)

## Replay the published figures

The paired plot subtracts `0_original` within each fold and seed. The figure
receipt records the plotted CSV, PNG and renderer hashes. [CSV-only procedure](AGENTS.md#replay-the-published-figures) ·
[Epoch traces](figure/data/learning-curve.csv) ·
[Held-out plot values](figure/data/metrics.csv) · [Paired plot values](figure/data/metrics-paired.csv) ·
[Figure receipt](log/figure-render.json)
