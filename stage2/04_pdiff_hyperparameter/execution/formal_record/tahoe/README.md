# Tahoe formal held-out record

Nine random-init PerturbDiff fits completed across folds 0–2 and seeds 17, 42 and 2026. Runs stopped at 7,000–9,000 optimizer steps and selected EMA checkpoints at 7,000–9,000 steps. The [frozen plan](plan.json) binds the Tahoe release, shared 18,164-gene CM axis, [validation-only pilot recipes](../../tahoe_pilot_record/README.md) and worker identities before held-out sampling.

The folds contain 18,896, 18,900, 18,731 held-out groups and 22,493, 21,850, 20,975 wells, respectively. The evaluator scored 9 PerturbDiff slots and six matched reference slots. Each reference was scored once per fold on the same frozen groups. The independent the original Git record matched every output byte.

The original retrieval score table and figures were removed. Active tenv
reporting uses `AUPRC`, `AUPRC_positive`, and `AUPRC_negative` through the
[shared evaluator](../../../../evaluation/README.md). Tahoe test retrieval is
outside the reporting scope. Training metadata and checkpoints remain intact.

| Scheduler phase | Completed jobs | Interrupted attempts | Observed span (min) | Peak host RSS (GiB) | Peak workstation GPU use (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| freeze | 1 | 0 | 0.4 | 1.105 | 1.773 |
| fit | 9 | 0 | 100.7 | 5.901 | 6.216 |
| predict | 9 | 2 | 406.6 | 4.177 | 3.559 |
| collect | 1 | 0 | 0.0 | 0.038 | 0.543 |
| evaluate | 1 | 0 | 16.1 | 15.081 | 0.549 |
| replay | 1 | 0 | 16.4 | 15.564 | 0.596 |
| audit | 1 | 0 | 0.1 | 0.030 | 0.555 |

The prediction span sums the observed intervals before and after the reboot.

The evaluation and replay jobs ran in scopes configured for 10 GiB of charged memory and 8 GiB of swap. Their `time -v` process RSS peaks were 15,813,064 and 16,319,664 KiB, respectively. The evaluator uses file-backed memory maps, and process RSS includes resident mapped pages. The scheduler recorded the scope caps and process RSS separately.

The power cut interrupted prediction for `fold-1-seed-17` and `fold-1-seed-42`. The scheduler preserved the original attempts and retried prediction from the sealed fits. The execution record retains both attempts for each recovered slot. The [recovery record](recovery.json) identifies the preserved partial output and archived telemetry.

Individual fits took 9.9–12.2 minutes. [Training provenance](training-provenance.json) records the scheduler timing, resource peaks, attempt logs and checksums. The slot directories preserve each original fit history and the verified prediction and checkpoint replay manifests. Raw telemetry and full artifacts remain local under `execution/output/tahoe-formal-ticket73/` and `execution/output/tahoe-formal-ticket73-control/`.

To reproduce the coverage audit from those local outputs, run `python stage2/04_pdiff_hyperparameter/execution/script/audit_tahoe_formal.py stage2/04_pdiff_hyperparameter/execution/output/tahoe-formal-ticket73 5148ddcd1f6df7da8b17265ccba0fad9cfe1831adcb8970d2fb0bf665354d556` from the repository root. The execution record pins the audit source and result hashes.
