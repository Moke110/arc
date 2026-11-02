### Evaluated population (test role)

- 951 rows / 347 conditions / 118 drugs over 5 folds; identical to #8's test population.
- Metric: direction AUPRC@200 (macro over drug×dose conditions, replicates aggregated within a condition first); auxiliary Δ RMSE / Pearson / Spearman over the same conditions.

### Headline test results by initialisation

| init | AUPRC@200 (macro over folds) | per-fold (mean of 3 seeds) | up | down | Δ RMSE | Pearson | Spearman |
|---|---|---|---|---|---|---|---|
| random | **0.0886** | 0.0910 / 0.0976 / 0.1033 / 0.0740 / 0.0770 | 0.0076 | 0.1696 | 0.3338 | 0.0680 | 0.0687 |
| tahoe | **0.0823** | 0.0852 / 0.0903 / 0.0952 / 0.0662 / 0.0745 | 0.0093 | 0.1552 | 0.4043 | 0.0773 | 0.0757 |

**Simple baselines on the identical test population (#8):**

| baseline | AUPRC@200 | up | down | Δ RMSE | Pearson | Spearman | role |
|---|---|---|---|---|---|---|---|
| `no_change` | 0.0110 | 0.0110 | 0.0110 | 0.1464 | n/a | n/a | primary |
| `delta_mean` | 0.1645 | 0.1572 | 0.1717 | 0.1419 | 0.2479 | 0.1989 | primary |
| `delta_mean_dose` | 0.1608 | 0.1508 | 0.1708 | 0.1424 | 0.2403 | 0.1886 | auxiliary |
| `ridge_fp` | 0.1512 | 0.1405 | 0.1618 | 0.1510 | 0.2175 | 0.1664 | auxiliary |

The primary bar is `delta_mean` = **0.1645** AUPRC@200.

### Acceptance question (a): does Tahoe pretraining beat random init?

- **Paired drug-level gain (Tahoe − random): -0.0063** AUPRC@200, 95% bootstrap CI [-0.0074, -0.0052], n=118 drugs, B=10000, seed=20261102.
- **Verdict (pre-registered rule): `negative`.**
- Seed stability (mean paired gain per seed): 20261021: -0.0028, 20261022: -0.0090, 20261023: -0.0070 (sd over seeds 0.0031).

| fold | mean paired gain over drugs |
|---|---|
| 0 | -0.0059 |
| 1 | -0.0072 |
| 2 | -0.0079 |
| 3 | -0.0077 |
| 4 | -0.0027 |

### Acceptance question (b): does either init beat the simple baselines?

- **Tahoe init − `delta_mean`: -0.0800** AUPRC@200 (95% CI [-0.1025, -0.0588], n=118 drugs) — verdict `negative`.
- **random init − `delta_mean`: -0.0738** AUPRC@200 (95% CI [-0.0965, -0.0521], n=118 drugs) — verdict `negative`.
- vs `no_change` (chance AUPRC = k/n_genes): Tahoe +0.0711 (CI [+0.0620, +0.0808]).

### Strata

| stratum | value | n conditions | n drugs | random AUPRC@200 | Tahoe AUPRC@200 | Tahoe−random | 95% CI |
|---|---|---|---|---|---|---|---|
| dose | 0.2 | 116 | 116 | 0.0824 | 0.0759 | -0.0062 | [-0.0074, -0.0051] |
| dose | 1 | 114 | 114 | 0.0868 | 0.0810 | -0.0064 | [-0.0076, -0.0053] |
| dose | 10 | 117 | 117 | 0.0967 | 0.0900 | -0.0062 | [-0.0074, -0.0051] |
| tahoe_seen | True | 45 | 15 | 0.0885 | 0.0816 | -0.0069 | [-0.0101, -0.0039] |
| tahoe_seen | False | 302 | 103 | 0.0886 | 0.0824 | -0.0062 | [-0.0074, -0.0049] |
| response_label | unknown | 347 | 118 | n/a | n/a | n/a | [n/a, n/a] |
| toxicity_label | unknown | 347 | 118 | n/a | n/a | n/a | [n/a, n/a] |

Low-response and toxicity labels are **unknown** on all 347 conditions: no per-well readout exists in this round and the labels are never inferred from transcriptomic magnitude (06.2 §6).
Tahoe-seen strata: 15 of the test drugs (45 conditions) have a drug key observed in Tahoe treated rows; 103 drugs (302 conditions) were never observed there. (354 drug keys are observed in Tahoe treated rows overall.)

### Actual budgets and failures

- Evaluation: 30/30 units completed, **0 failed**; wall 208.4 s total, per-unit mean 6.9 s; host peak 897 MiB; GPU peak 1593 MiB; repeated inference byte-identical (predict.py runs it twice).
- random training budget of the evaluated checkpoints: epochs 27/67.7/100 (min/mean/max), capped=3, stop reasons `{'early_stop_no_auprc_improvement_20': 12, 'max_epochs': 3}`.
- tahoe training budget of the evaluated checkpoints: epochs 30/44.9/100 (min/mean/max), capped=1, stop reasons `{'early_stop_no_auprc_improvement_20': 14, 'max_epochs': 1}`.

### Input and model identity

- units manifest: `stage2/06_prnet_transfer/06.9_evaluation/output/record/units.json` sha256 `62b79b622622868f…` (created 2026-11-02T11:34:51Z, 30/30 units scored).
- frozen data identity: `a57b0e35b1126d0b88da361a43ff27479f42e005a5291cd7641d3fa52cb7b376` (manifest sha256 `3050d149198cb3d2…`).
- every checkpoint sha256 and config digest re-verified against #6/#7's `selected_checkpoints.csv` before scoring; model cards in `output/model_cards/`.

### Recorded deviations

- `variance_head_guard_relaxed` (stage2/06_prnet_transfer/06.4_pipeline/execution/ptlib/model.py::check_architecture, fix `813c0bb`): The all-zero architecture probe required the decoder variance to be strictly positive. Four selected Tahoe-init checkpoints underflow softplus to exactly 0 on that off-distribution input; the guard now rejects only non-finite or negative variance and accepts an exact zero. Units: e9_f0_s20261023_lr1e-4_tahoe, e9_f1_s20261023_lr1e-4_tahoe, e9_f4_s20261022_lr1e-4_tahoe, e9_f4_s20261023_lr1e-4_tahoe. Effect on results: none on any predicted value: guard-only change, no model, loss, checkpoint or selection changed; the 26 already-scored units are byte-identical (eval manifest + unit record hashes) and only the 4 failed units were re-scored; 4/4 re-scored successfully; final 30/30 units completed.

### Scope, negative-result clause and decision rule

- **Scope.** This round predicts general 24 h iPSC-CM transcriptomic responses to drugs never measured in the temp_env CM training set, on the single fixed training-consensus background c*. It evaluates the 30 checkpoints selected before the test split was scored. It is not a hypoxia/reoxygenation efficacy claim, not a ScreenSeq or DToxS result, not a dose-extrapolation claim (the three measured temp_env doses 0.2/1/10 µM are reported as strata), and not a donor- or batch-specific claim. Low-response and toxicity labels are unknown in this round and are reported as unknown.
- **Negative-result clause.** A negative result is a valid first-round outcome. If the paired Tahoe-vs-random gain or the advantage over the simple baselines is not positive (95% bootstrap CI over drugs excluding zero), the study reports that no transfer gain was demonstrated and that PRnet on this frozen input should not be described as a drug-screening tool. No tuning, re-selection or seed/fold/LR change is performed after seeing the test result; the 30 units are the frozen #6/#7 selection.
- **Decision rule.** Pre-registered reading: the primary evidence is the drug-level paired mean difference in direction AUPRC@200 with a 95% percentile bootstrap CI resampling drugs (B=10000, seed=20261102). 'positive' = point estimate > 0 and CI low > 0; 'negative' = point estimate <= 0; otherwise 'inconclusive'. Both acceptance questions (vs random init; vs the simple baselines) are reported regardless of sign, and all three seeds and five folds are shown — no best-seed selection.

### Reproduce

```bash
cd stage2/06_prnet_transfer/06.9_evaluation/execution
make units            # verify the 30 frozen selections (checkpoint sha256/config/identity)
make score-test       # predict + evaluate all 30 units (test phase, explicit confirmation)
make record           # tables + summary + gates  -> output/record
make cards            # 30 model cards            -> output/model_cards
make report           # regenerate this section + figures
make check            # re-derive the record byte-for-byte
```
