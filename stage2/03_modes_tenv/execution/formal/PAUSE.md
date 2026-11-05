# Formal tenv execution pause — 2026-11-06

Paused at the user's explicit request after the active computation finished
naturally. **34/63 formal runs completed, zero failures, 29 not dispatched.**
No formal test scoring or prediction-barrier sealing has occurred. All model
admission tickets are complete; this pause does not reopen admission.

## Durable state

- Integration branch: `dev`; code/registry/preflight commit before pause: `9185aa8`.
  The subsequent pause commit contains handoff documents and compact evidence.
- Completed: `no_change` 9, `delta_mean` 9, PRnet 9, CPA 7.
- Remaining: CPA fold2/seed42 and seed2026; chemCPA 9; CMonge 9; PerturbDiff 9.
- Exact registry, commands, 63 expected identities and unchanged serial queue:
  `registry.json`, `commands.json`, `expected_identities.json`, `queue.json` here.
- `pause_record.json` records verification of all frozen identities and completed
  artifacts. `pause_record/completed/` holds compact run/checkpoint/scheduler
  copies. Full prediction arrays and checkpoint binaries remain local-only under
  `stage2/03_modes_tenv/output/formal/`; metadata copies are not binary uploads.
- The last worker, CPA fold2/seed17, naturally exited 0. Native training 216.393 s,
  inference 0.733 s; wrapper wall 237.09 s, peak RSS 1,871,652 KiB. Its completed artifacts
  were verified before scheduler state was finalized from checkpoint times and
  the `/usr/bin/time -v` log, following the scheduler's success-accounting policy.
- To prevent new dispatch without interrupting training, only dispatcher 736317
  was suspended. It was removed after worker 747319 naturally exited 0; no model
  worker was killed. Dispatcher sampling was suspended during this requested
  drain; manual process/resource checks and the final drain-guard record are
  retained. The worker's scope cap and adapter time limit remained active.
- No experiment scheduler, worker, adapter or drain monitor remains active.
  No merge to main or push was performed.

## Resume

1. Read the latest #40/#46/#47 pause comments, this file, and `pause_record.json`.
   Confirm `git status --short --branch` is clean on `dev`. Read the execution
   environment rules and inspect actual host memory/GPU/process headroom.
2. Preserve the frozen input release, source files, environments and configs.
   Identity checks include installed package metadata; CMonge's CUDA overlay
   files and explicit reused-library roots are verified live. Missing or stale
   binaries must not be silently replaced or identities manually rewritten.
3. Run the **same original queue**, on the actual host rather than a restricted
   PID namespace. The scheduler validates/skips the 34 completed runs and starts
   the 29 untouched jobs. Do not delete scheduling state, reset retry budgets,
   copy admission pilots into formal runs, or create a second scheduler.

```bash
/home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/schedule.py run \
  stage2/03_modes_tenv/execution/formal/queue.json
```

The queue has `max_active_jobs: 1`, at most four hours of cumulative training per run,
separate inference limits, one environment repair/retry and one OOM batch retry.
The scheduler provides five-second telemetry and 10 GiB/8 GiB scope caps. Expected long-work
projections remain approximately 18 h for nine native CMonge runs and 30 h for nine
native PerturbDiff runs including inference; projections are not guarantees.

4. Keep all test scoring withheld until all 63 outcomes are completed or genuine
   terminal failures/timeouts. Then seal using the exact frozen model set and
   command mapping:

```bash
/home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/benchmark.py seal \
  --input stage2/03_modes_tenv/output/input \
  --output stage2/03_modes_tenv/output/formal \
  --models no_change delta_mean prnet cpa chemcpa cmonge perturbdiff \
  --commands stage2/03_modes_tenv/execution/formal/commands.json
/home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/benchmark.py evaluate \
  --input stage2/03_modes_tenv/output/input \
  --output stage2/03_modes_tenv/output/formal --workers 2
/home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/benchmark.py replay \
  --output stage2/03_modes_tenv/output/formal
```

5. Complete #47 tables, source arrays and three PNG/SVG figures; retain undefined
   counts, deterministic reference copies, post-hoc best-test labels and genuine
   incomplete outcomes. Run final offline health gates under resource limits
   after heavy training ends. Perform independent Standards/Spec code-review
   against `f022b1b`, resolve material findings, commit on `dev`, then close
   #46/#47/#40 with final evidence. Do not push `dev`.

## Checks already passed

Five predictors passed current full 18,164-gene runtime/resource admission:
PRnet, CPA, chemCPA, CMonge and PerturbDiff. The final registry retains 11
excluded and 4 unresolved candidates with evidence. Native objectives and
interfaces, actual installed environments, selected config hashes and input
separation are recorded in each model's fidelity/admission files.

Before dispatch: 33/33 offline provenance verifiers pass, path/tracked guards
pass, pipeline 4 / metrics 3 / environment-identity 3 / scheduler 6 tests pass. Independent
core preflight findings were fixed and re-audited; the separate CMonge/
PerturbDiff static admission audit found no actionable P1/P2 issue. These checks
are not the pending final post-matrix review or scientific evaluation.

Admission tickets #45/#48/#49/#50/#51 are closed. Implementation #41–#44 were
already closed. **#40, #46 and #47 stay open for the requested continuation.**
