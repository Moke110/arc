# Step 04 — PerturbDiff training configuration

Step 04 studies standalone PerturbDiff training on Tahoe and tenv on the
18,164-gene CM axis. [ADR 0046](../../docs/adr/0046-perturbdiff-standalone-training.md)
defines the validation, checkpoint and reporting contract.

## Queue operation

The dedicated scheduling subagent runs the queue from the repository. The main
agent relays user controls and can read status during execution.

```bash
scheduler=stage2/04_pdiff_hyperparameter/execution/script/schedule.py
python "$scheduler" preflight "$queue"
python "$scheduler" run "$queue"
python "$scheduler" status "$queue"
python "$scheduler" pause "$queue"
python "$scheduler" resume "$queue"
python "$scheduler" retry "$queue" JOB_ID
```

Use a frozen JSON queue with `record_dir` and a `jobs` array. Each job has a unique
`id`, an argument-list `command`, a `kind` of `cpu`, `gpu` or `training`, and an
output-space estimate `disk_gib` (default 1). Set `profile.host_peak_gib` and
`profile.gpu_peak_gib` from a measured run. A first profile without measured
peaks starts alone and reserves the worker's full RAM cap and at least 4 GiB
of GPU memory for GPU work, before headroom. `cwd` defaults to the scheduler's
repository root. Relative paths for `record_dir` and `cwd` resolve against the
queue file. The record directory must share the filesystem of the job's outputs
for the free-disk gate to apply to those outputs. `timeout_seconds` defaults to
14,400. Use `kind: training` for exclusive long fits. `memory_gib` defaults to 10
and permits 18 for a formal exclusive fit. Swap is capped at 8 GiB.

Preflight is read-only. A run holds the Step 03 workstation lock, samples live
host RAM, disk and GPU every two seconds, and starts at most one worker per
sample. Admission requires 12 GiB host RAM to remain after the predicted job
peak plus 20% headroom, at most two heavy workers, approximately 28 GiB combined
peak, and GPU and disk headroom. The Step 03 process scan detects GPU compute
clients and large Python, R or Julia processes. Visible external heavy work
blocks dispatch because its future peak is unknown. Three consecutive unsafe
samples stop an owned worker. The GPU admission gate uses utilization at or
above 90% for new GPU jobs. In-flight GPU pressure uses free VRAM because a
healthy training worker can keep the GPU busy.

Pause lets active workers finish and holds the queued jobs. After the queue
drains, resume clears pause and the dedicated scheduler can invoke run again.
The scheduler preserves interrupted attempts and requires explicit retry.
When a previous scheduler died, run or retry verifies that its named scope and
process have exited before changing the outcome. A live orphan blocks recovery
and must drain before retry. Keep one immutable queue per record directory.

`status.json` records worker identity, attempt outcomes and peaks.
`telemetry.jsonl` records live decisions' input samples. Worker logs retain
`/usr/bin/time -v` peak RSS on normal exit. RSS sampled across the owned process
tree and GPU use sampled across the workstation are recorded separately from
the resource caps. A completed scheduler worker means its command exited
successfully. Its experiment artifacts still go through the corresponding
artifact verifier.

## Full-axis tenv smoke

The bounded Step 04 smoke uses fold 0, seed 17 and the frozen Step 03 tenv
release. The preparation command pins the release to the Step 03 verification
record and verifies the ordered gene axis, split and predictor files by hash.
The resulting request binds the consumer manifest and allowlisted files.
Predictors receive held-out covariates for the frozen drug and dose vocabulary.
The smoke uses a validation group for its prediction probe. Test-treated
expression and labels remain in the evaluator directory.

From the repository root, prepare the smoke request in a fresh scratch
directory under the step. Give the GPU commands to the dedicated scheduler as
separate queue jobs.

```bash
step=stage2/04_pdiff_hyperparameter
cli="$step/execution/script/train.py"
pdiff=stage2/_envs/perturbdiff/bin/python
smoke="$step/execution/output/tenv-smoke-replay"
mkdir -p "$smoke"
"$pdiff" "$cli" prepare-tenv stage2/03_modes_tenv/output/input "$smoke/request.json" --fold 0
"$pdiff" "$cli" verify-request "$smoke/request.json"

# First scheduler job.
"$pdiff" "$cli" fit "$smoke/request.json" "$step/execution/smoke_record/recipe.json" "$smoke/run"

# Follow-up scheduler job after the fit exits successfully.
"$pdiff" "$cli" reload "$smoke/request.json" "$smoke/run" "$smoke/replay.json"
```

`fit` requires a fresh run directory. It records the actual fixed model
recipe and validation batching in `run.json`. Every validation event ranks the
current EMA state by strict minimum of the native monitor. Patience measures
optimizer steps since the most recent improvement strictly greater than
`min_delta`. A budget or signal interrupt leaves an incomplete run record
when the worker can exit gracefully. `best.pt` stores the selected EMA
prediction state. The checkpoint timing includes producing the validation
prediction probe and writing and hashing the checkpoint. `reload` verifies
the checkpoint and model state hashes, recomputes the selected validation
loss and runs the same native predictor with identical probe noise.

The smoke kept all 138 validation groups, with 16 full batches of eight
groups and one final batch of two. It validated and saved at steps 6 and 12,
ending at its planned 12-step cap. The selected validation loss was
14.8742036819458. The isolated reload reproduced the loss and the
3 × 18,164 prediction exactly. The worker's peak RSS was 3,049,412 KiB.
See [the run record](execution/smoke_record/run.json) and
[scheduler evidence](execution/smoke_record/execution.json). The original
checkpoint and raw scheduler telemetry are local under
`execution/output/tenv-smoke-ticket67/`. The smoke is outside the formal
matrix.
