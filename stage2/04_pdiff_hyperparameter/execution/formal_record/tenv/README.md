# Tenv formal held-out record

Nine random-init PerturbDiff fits completed across folds 0–2 and seeds 17,
42 and 2026. All stopped by the selected validation policy at 5,500–8,000
optimizer steps. The selected EMA checkpoints occur at 5,000–7,500 steps.
The [frozen plan](plan.json) binds the Step 03 tenv release, shared
18,164-gene CM axis, the [pilot recipes](../../pilot_record/tenv/README.md),
predictor requests and launch identities before any held-out sampling.

The evaluation reuses the held-out drug memberships already scored in
Step 03. The three folds contain 462, 459 and 459 held-out groups across
1,401, 1,407 and 1,398 wells, respectively. The two references are scored
once per fold on the same groups as PerturbDiff.

The original retrieval score table and figures were removed. Active tenv
reporting uses `AUPRC`, `AUPRC_positive`, and `AUPRC_negative` through the
[shared evaluator](../../../../evaluation/README.md). Tahoe test retrieval is
outside the reporting scope. Training metadata and checkpoints remain intact.

| Scheduler phase | Completed jobs | Elapsed (min) | Peak host RSS (GiB) | Peak workstation GPU use (GiB) |
| --- | ---: | ---: | ---: | ---: |
| Fit | 9 | 63.7 | 3.037 | 6.205 |
| Held-out prediction | 9 | 14.3 | 2.478 | 2.831 |
| Evaluation | 1 | 0.5 | 1.007 | 1.832 |
| Checkpoint release | 1 | 0.1 | 0.038 | 1.738 |

Elapsed time is the phase's monotonic telemetry span. Each job completed on
its first attempt. [Training provenance](training-provenance.json) records per-slot
timings, worker and native resource peaks, scheduler attempts and hashes
for status, logs and raw telemetry. The original records, predictions and
group-level scores are local at `execution/output/tenv-formal-ticket72/`.

All nine downstream checkpoints are local at `output/checkpoints/tenv/`.
Use the [checkpoint index](checkpoints.json) to select the exact fold and
seed. Every inference worker reloaded the selected EMA, reproduced the
validation loss exactly and matched the saved validation prediction with
zero maximum absolute error before opening any held-out batch. Each slot
directory in this record contains the original compact run history,
prediction manifest and replay evidence. The scheduler's publication job
verified the checkpoint hashes and hard-linked each published checkpoint
to its fitted artifact. The [Step 04 README](../../../README.md#select-a-downstream-tenv-checkpoint)
gives the downstream replay command.
