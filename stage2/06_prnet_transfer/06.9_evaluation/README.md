# 06.9_evaluation — 30-unit test evaluation and the first-round report (ticket #9)

Status: **harness built, offline-verified, and dry-run on the real validation role; the 30
test units are pending the Lead's freeze signal.** Ticket #9 is the only place the test
split of `06.3_input` is ever scored for the first-round PRnet study. The 30 units are
exactly the checkpoints **selected before test** by
[`06.6_cm_random`](../06.6_cm_random/README.md) (15 random-init) and
[`06.7_cm_tahoe`](../06.7_cm_tahoe/README.md) (15 Tahoe-init); `06.8_baselines` supplies the
paired `delta_mean` / `no_change` reference on the identical test population.

- decision: [ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md)
- design: [`../README.md`](../README.md) · spec `.scratch/stage2-prnet-formal-transfer/spec.md`
- frozen input: [`06.3_input`](../06.3_input/README.md)
- evaluator / predictor / model cards: [`06.4_pipeline`](../06.4_pipeline/README.md) — **reused, never forked**
- baselines: [`06.8_baselines`](../06.8_baselines/README.md)

## 1. What is delivered

| deliverable | where (after the test run) |
|---|---|
| the 30-unit manifest with **checkpoint sha256 + config digest + data identity**, fully re-verified before scoring | `output/record/units.json`, `output/record/units_check.json` |
| every unit's test prediction/evaluation (frozen evaluator output) | `output/units/<unit_id>/{pred,eval}/` (local, reproducible) |
| per-unit, per-fold×seed and per drug×dose results | `output/record/unit_results.csv`, `fold_seed_table.csv`, `per_drug_dose.csv` |
| paired Tahoe-vs-random gain **by drug** with bootstrap uncertainty and seed stability | `output/record/paired_by_drug.csv`, `drug_gains.csv` |
| dose / Tahoe-seen / unknown low-response / unknown toxicity strata | `output/record/strata.csv` |
| headline summary, budgets/failures, identity, scope + negative-result clause | `output/record/summary.json`, `manifest.json`, `record_check.json` |
| 30 model cards (the #4 writer, pointed at each run + its evaluation) | `output/model_cards/*.md` + `index.json` |
| figures | `output/figures/*.png` |
| this report's results block | generated section below (§7), from the tracked tables |
| reproduce entrypoint | `execution/run.sh` / `execution/Makefile` (§5) |

Everything scientific is derivable from the frozen predictions: the tables are rebuilt by
`write_record.py` from `per_condition.csv` / `aggregate.json` / `deltas.npz` that the frozen
evaluator wrote, and `make check` re-derives the record byte-for-byte.

## 2. The 30 units and their verification

`5 folds × 3 seeds × 2 initialisations`, one unit per (fold, seed, init) at that branch's
**selected common LR** — the selection rule #6/#7 locked (three-seed validation-AUPRC mean,
lower LR on a tie). No seed, fold or candidate is dropped or cherry-picked.

`build_units.py` fails closed unless all of the following hold, and records each as a gate:

1. **Selection records.** One `selected_checkpoints.csv` per branch (15 rows each), exactly
   5 folds × the 3 frozen seeds, one common LR per (fold, init), one data identity.
2. **Checkpoint bytes.** Every selected checkpoint resolves (searching
   `$ARC_CHECKPOINT_ROOTS`, the worktree, then the main checkout — the weights are
   gitignored local files), and its **sha256 + size** match the record.
3. **Config digest chain.** Per candidate run:
   `sha256(config_fingerprint(config.json)) == config.json.config_sha256 ==
   run manifest.config_sha256 == selected record.config_sha256`, and the run's fold/seed/
   LR/structure/metric-k/train-set agree with the record and the locked config.
4. **Data identity.** `content_sha256({manifest_sha256, locked_config_sha256, artifacts})`
   equals the identity every selected row and candidate manifest records; all 30 units share
   it. (Verified on the real record: `a57b0e35b1126d0b88da361a43ff27479f42e005a5291cd7641d3fa52cb7b376`.)
5. **No test contact during selection.** No candidate run directory holds a
   `prediction/`, `evaluation/` or `test/` artifact and no finetune manifest records a role.
6. **Init provenance.** Random units: `init.mode == "random"`. Tahoe units:
   `init.mode == "pretrained"` with the #5 warm-start checkpoint sha256 recorded; if the
   warm-start file is still on disk its sha256 is re-verified (a missing gitignored weight is
   reported as a warning, never as "verified").

Before scoring, `score.py` additionally checks the produced prediction/evaluation manifests:
the evaluator's `predictions.checkpoint_sha256` equals the unit's verified checkpoint, the
evaluator role is the requested role, the prediction array scored is the one predicted, the
repeated inference is byte-identical, and the scored row/condition counts equal #8's recorded
population for that fold (192/188/191/190/190 test rows and 68/67/71/71/70 conditions).

## 3. Scoring (frozen evaluator only)

Per unit, two ticket-#4 CLIs are invoked as subprocesses, never re-implemented:

```bash
predict.py  --stage finetune --fold F --role {valid,test} --checkpoint <verified> --repeats 2
evaluate.py --stage finetune --fold F --role {valid,test} --pred-dir <pred>
```

- metric: **direction AUPRC@200**, macro over `drug_key × dose` conditions after
  within-condition replicate aggregation; auxiliary Δ RMSE / Pearson / Spearman.
- truth Δ: experimental treated log1p − per-plate matched DMSO mean, read **only** by the
  evaluator, never a model input (validation/test inference gets the fold's fixed `c*`).
- test is scored **once**, after candidate selection, and never feeds back into LR, epoch,
  seed or capacity. `--phase test` refuses to run without `--confirm-frozen-selection`; the
  Makefile target is additionally scoped with `systemd-run` and named for that gate.
- one unit at a time; per-unit wall time, host peak RSS and whole-GPU peak are recorded;
  repeated inference byte-identical is a hard requirement, not a warning.
- estimated cost: the full test sweep is ~30 × 7–10 s ≈ **4–6 min**, < 1 GiB host, < 3 GiB
  GPU — far below the 8.5 GiB full Tahoe-valid evaluation, and it never materialises an
  observation × 18,164 dense matrix.

## 4. Aggregation and uncertainty

- **Per fold/seed:** every unit is reported; fold×init tables carry the three-seed mean/sd.
- **Per drug×dose:** the evaluator's per-condition table for all 30 units, joined (for the
  identical population) to #8's `no_change`, `delta_mean`, `delta_mean_dose`, `ridge_fp`
  AUPRC@200 per condition.
- **Paired Tahoe-vs-random gain by drug.** Within `(fold, seed, drug)` the arm score is the
  equal-weight mean over that drug's dose conditions; the paired difference is
  `tahoe − random`; the per-drug gain averages over the three seeds. Doses, cell lines and
  candidate runs are **never** treated as independent drugs. Uncertainty is a percentile
  bootstrap resampling **drugs** (B = 10000, seed 20261102, 95% interval); the per-seed and
  per-fold gains and the seed sd are reported alongside.
- **Baseline comparison.** The same pairing against #8's primary `delta_mean` bar and the
  `no_change` constant, on the identical 951-row / 347-condition test population.
- **Strata.** Measured doses (0.2/1/10 µM); Tahoe-seen drug keys (observed in Tahoe treated
  rows — a structure-seen flag, not a CM label). Low-response and toxicity labels are
  **unknown on every condition this round**: 06.2 §6 establishes there is no per-well
  readout and forbids inferring them from transcriptomic magnitude. They are reported as an
  explicit `unknown` stratum with that pre-registered criterion, never dropped or guessed.
- **Pre-registered reading:** `positive` = paired mean gain > 0 and 95% CI low > 0;
  `negative` = point estimate ≤ 0; otherwise `inconclusive`. Both acceptance questions
  (vs random init; vs the simple baselines) are reported regardless of sign.

## 5. Reproduce

```bash
cd stage2/06_prnet_transfer/06.9_evaluation/execution

make test                     # 47 offline unittest tests (no data, no GPU)
make units                    # THE 30-unit manifest + verification (fails closed until #7
                              # record exists; --allow-incomplete is the harness dry run)
make plan                     # print the scoring plan; reads no test row
make score-valid              # harness smoke: validation role only (no test data)
make score-test-scoped        # THE single test touch (systemd-run 8 GiB cap; explicit gate)
make record                   # tables + summary + gates -> output/record
make cards                    # 30 model cards            -> output/model_cards
make report                   # figures + results block   -> README §7
make check                    # re-derive the record, byte-compare
```

`score-test-scoped` is the only target that sets the frozen-selection confirmation; run it
only after the Lead confirms both matrices are final. `run.sh` exposes the same commands.

Defaults: `DATA_DIR=/home/chang/arc/stage2/06_prnet_transfer/06.3_input/output`,
`BASELINES_DIR=…/06.8_baselines/output`, `TABLE_CACHE=…/06.4_pipeline/output/input_cache`
(the pinned env has no parquet engine; build it once with `make tables` if absent), and the
records default to `06.6_cm_random`/`06.7_cm_tahoe` under the current tree. Override
`T09_DATA_DIR` / `T09_OUT` / `ARC_CHECKPOINT_ROOTS` as needed.

## 6. Harness validation (no test science)

- **Offline suite:** 47 tests over synthetic units/eval directories — identity-chain
  tampering (sha256, size, config digest, data identity), missing/extra/duplicate units,
  LR-rule violations, test-artifact detection, macro recomputation, per-drug pairing
  (doses not counted as drugs), bootstrap determinism, strata, baseline join, dry-run
  labelling, byte-identical record re-derivation, report splicing.
- **Real validation dry run** (validation role only, `output/dry_run/`, never tracked):
  fold-0 seed-20261021 random-init unit scored end-to-end through the real CLIs —
  188 valid rows / 67 conditions, 6.9 s, **887 MiB** host peak, repeated inference
  byte-identical, and the evaluator's macro AUPRC@200 (0.0968342) reproduced the
  checkpoint's training-time best validation AUPRC to 1e-9. The aggregation, gates and
  report rendered a correctly labelled dry-run record. This exercises the full harness on
  the real frozen input without touching the test split.

- **Recorded deviation — variance-head guard (commit `813c0bb`).** The frozen #4
  `check_architecture` probe (all-zero batch) required the decoder variance to be strictly
  positive. Four selected Tahoe-init checkpoints (`e9_f0_s20261023_lr1e-4_tahoe`,
  `e9_f1_s20261023_lr1e-4_tahoe`, `e9_f4_s20261022_lr1e-4_tahoe`,
  `e9_f4_s20261023_lr1e-4_tahoe`) underflow `softplus` to exactly 0 for 464/300/2/482 of the
  36,328 probe values, because float32 underflows at very negative logits. Their weights are
  finite, on real test inputs `var.min` is 2.19e-12 / 3.80e-11 / 8.79e-13 / 1.16e-8 with
  **zero** zeros, and the evaluated metric uses the **predicted mean only**. The guard was
  relaxed (reject non-finite or negative variance only) with six regression tests in
  `06.4_pipeline/execution/tests/test_model_guard.py`; no model, loss or predicted value
  changed. The 26 already-scored units are byte-identical before/after (eval manifest + unit
  record hashes) and only those 4 units were re-scored (4/4 recovered, final 30/30). The
  deviation and its evidence are recorded in `output/record/summary.json` and
  `manifest.json` (`deviations[0]`).

## 7. First-round results

<!-- BEGIN RESULTS (generated by 06.9_evaluation/execution/report.py; do not edit) -->
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

- units manifest: `/home/chang/arc-wt/t09-evaluation/stage2/06_prnet_transfer/06.9_evaluation/output/record/units.json` sha256 `b4174c8407489f70…` (created 2026-11-02T11:34:51Z, 30/30 units scored).
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
<!-- END RESULTS -->

## 8. Boundaries and notes

- **Negative result clause.** A negative result is a valid first-round outcome. If the paired
  Tahoe-vs-random gain or the advantage over the simple baselines is not positive (95%
  bootstrap CI over drugs excluding zero), the study reports that no transfer gain was
  demonstrated and that PRnet on this frozen input is not described as a drug-screening
  tool. No tuning, re-selection or seed/fold/LR change happens after seeing the test result.
- The test split is scored exactly once. `record_check.json` stores the gates that prove the
  scoring followed the frozen selection (manifest created before scoring, one evaluation
  per unit, test role recorded, identical population).
- Checkpoints and per-unit prediction/evaluation arrays are gitignored local output; the
  tracked record is `output/record/` (force-added), `output/model_cards/` and
  `output/figures/`. The record is re-derivable with `make check`.
- One recorded deviation exists (variance-head guard, commit `813c0bb`, §6): guard-only,
  no predicted value changed, four units re-scored after the fix; the other 26 unit
  evaluations are byte-identical. No selection, LR, epoch or seed changed.
- GitHub issue #9 stays open until the Lead's integration; this substep does not merge,
  push or close it.
