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

A paused GPU prediction queue can resume with `--admission-profiles` on
preflight and run. The JSON sidecar binds the original queue hash, a hashed
local snapshot of a completed comparable worker, and conservative measured
host and GPU peaks for named, otherwise unprofiled GPU jobs. It changes only
resource decisions. Commands, scope limits and job identities come from the
original queue. The scheduler records the sidecar and each affected attempt's
profile. Combined GPU reservations must fit before an overlap is admitted.

`status.json` records worker identity, attempt outcomes and peaks.
`telemetry.jsonl` records live decisions' input samples. Use `observed_monotonic`
for elapsed intervals within the same workstation boot. Wall timestamps remain
for correlation with the system journal. Worker logs retain
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

## Full-axis Tahoe smoke

The Tahoe release projects the corpus v2 native `log1p_cp10k` layer by Ensembl
identity onto the ordered CM axis. The historical FAM153B slot is zero.
Each treated profile's target subtracts the mean of the declared `DMSO_TF`
profiles for its plate and cell line. The release records all eligible case
IDs, matched control indices and stratum membership hashes, case-level
exclusions, source hashes and the split assignments. Each outer fold holds
out entire normalized parent-connectivity blocks. The training and
evaluation group is drug × dose × cell line.

From the repository root, give each command to the dedicated scheduler in
order. Use a fresh release and run directory.

```bash
step=stage2/04_pdiff_hyperparameter
cli="$step/execution/script/train.py"
pdiff=stage2/_envs/perturbdiff/bin/python
smoke="$step/execution/output/tahoe-smoke-replay"
release="$smoke/release"
"$pdiff" "$cli" build-tahoe "$release"
"$pdiff" "$cli" verify-tahoe "$release"
"$pdiff" "$cli" prepare-tahoe "$release" "$smoke/request.json" --fold 0
"$pdiff" "$cli" fit "$smoke/request.json" "$step/execution/tahoe_smoke_record/recipe.json" "$smoke/run"
"$pdiff" "$cli" reload "$smoke/request.json" "$smoke/run" "$smoke/replay.json"
```

Reserve 25 GiB for the build's bounded CSR projection and temporary memmap.
The verifier checks the pinned source hashes, ordered genes, source-derived
control pairing, split memberships, file inventory and all-zero slot.
Request preparation uses the predictor file allowlist. Per-fold test group
means and their directional labels live in `fold_N/evaluator/`.

The smoke covers fold 0 and seed 17. Its full validation order is first
appearance in the frozen validation case CSV. Eight groups per batch and
three profiles per group reproduce the fixed native recipe. Sampling is
without replacement where possible, and the final short batch contributes
one scalar to the native monitor.

The real release kept 65,318 of 65,618 treated profiles across 377 drugs.
Missing structures excluded 150 profiles each for Sacubitril/Valsartan and
Verteporfin. All 1,400 DMSO_TF profiles contributed to 700 control strata.
The build took 128 seconds and peaked at 7,848,152 KiB RSS. See the
[release manifest](execution/tahoe_smoke_record/input.json) for coverage by
fold and source and file identities.

The fit validated at steps 6 and 12 on 712 full batches and one final batch
of three groups. Step 12 had the lower native loss, 15.660353660583496.
The isolated GPU reload reproduced that loss and the 1 × 18,164 validation
prediction exactly. The fit worker peaked at 3,582,104 KiB RSS.
The [run record](execution/tahoe_smoke_record/run.json) retains the loss
history and validation group order hash, the
[reload record](execution/tahoe_smoke_record/replay.json) binds the prediction
to the checkpoint, and the
[scheduler record](execution/tahoe_smoke_record/execution.json) identifies
the original telemetry. The original checkpoint and full group IDs live
under `execution/output/tahoe-smoke-ticket68/`. This smoke is outside the
formal matrix.

## Tenv validation pilot

The validation-only tenv pilot chose a 500-step interval in all three folds.
Validation and checkpoint overhead was 10.81–11.34% of optimizer time at 500
steps and 21.56–22.98% at 250 steps. The frozen per-fold recipes use a
4,000-step minimum, 1,000-step patience and fold-specific improvement
thresholds. See [the pilot outcome](execution/pilot_record/tenv/README.md).

## Tahoe validation pilot

The Tahoe pilot chose a 500-step interval and a 4,000-step minimum in all
three folds. Validation and checkpoint overhead was 47.71–50.76% of
optimizer time at 500 steps and 94.83–97.83% at 250 steps. The frozen
recipes use 2,000-step patience in fold 0 and 1,000-step patience in folds
1 and 2, with fold-specific improvement thresholds. See [the Tahoe pilot outcome](execution/tahoe_pilot_record/README.md).

## Formal fits and held-out scoring

The [tenv formal record](execution/formal_record/tenv/README.md) reports
the completed matrix, matched references and checkpoint handoff.

Use one formal plan per dataset. Freeze it from the verified compact pilot
record and the matching predictor requests in fold order. The freeze command
prints the plan's SHA-256 and creates separate nine-slot queues for fits and
predictions. Keep that hash with the scheduler records. The plan pins the
requests, raw and normalized recipes, pilot selection, execution code and
runtime identities. Use a stable checkout for the duration of the queues.

```bash
step=stage2/04_pdiff_hyperparameter
execution="$step/execution"
cli="$execution/script/formal.py"
pdiff=stage2/_envs/perturbdiff/bin/python
formal="$execution/output/tenv-formal-replay"
requests="$execution/output/tenv-pilot-ticket70/requests"
"$pdiff" "$cli" freeze "$execution/pilot_record/tenv" "$formal" \
  "$requests/fold-0.json" "$requests/fold-1.json" "$requests/fold-2.json"
```

The dedicated scheduler runs `fit-queue.json` to its terminal outcomes,
then runs `predict-queue.json`. If needed, use the scheduler's retry command
on the relevant queue and slot ID, such as `fold-0-seed-17`. An interrupted
fit restarts from random initialization in a fresh attempt directory. A
prediction retry can reuse the sealed fit. A completed prediction is published
only after the worker reproduces the native validation loss and prediction
probe exactly, samples all held-out groups, verifies case coverage and
flushes the float32 prediction. Inference keeps each group together and
uses at most 24 profiles across groups in a batch. A group larger than that
runs as a single group. Failed attempts retain their evidence.

After both phases reach terminal outcomes, collect the slot inventory
and run the evaluator with the `arc` interpreter. Give full-size scoring
and checkpoint release to the dedicated scheduler. Use the SHA-256 printed
by the freeze command for `plan_sha`.

```bash
python "$execution/script/formal_outcomes.py" "$formal/plan.json" "$plan_sha" "$formal/outcomes.json"
python "$execution/script/evaluate.py" "$formal/plan.json" "$formal/outcomes.json" "$formal/evaluation" \
  --plan-sha256 "$plan_sha"
python "$execution/script/formal_outcomes.py" "$formal/plan.json" "$plan_sha" "$formal/outcomes.json" \
  --publish-checkpoints "$formal/evaluation" "$step/output/checkpoints/tenv"
```

The [evaluation contract](execution/evaluation_contract.md) describes the
score tables and provenance. The standalone tenv evaluation uses the
held-out memberships previously scored in 03. The checkpoint release
requires all nine tenv slots to be scored and uses hard links on the same
filesystem for the large `best.pt` files.

### Select a downstream tenv checkpoint

Read `output/checkpoints/tenv/checkpoints.json` and select the exact
`identity.fold` and `identity.seed` for the intended split and initialization.
Each entry names a relative directory and records the frozen request,
recipe, checkpoint and model state hashes. The directory contains the
checkpoint, request, recipe, training record and exact EMA replay evidence.
To verify a selected checkpoint against its frozen validation input, give
this command to the dedicated scheduler with `slot` set to that directory.

```bash
"$pdiff" "$execution/script/train.py" reload "$slot/request.json" "$slot" "$fresh_replay_record"
```

The request retains the absolute path to the local tenv release and verifies
its allowlisted files before loading the model. The replay validates the
selected loss and prediction probe using the saved EMA weights.
