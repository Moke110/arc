# 03_modes_tenv — tenv shared-Morgan comparison

This step implements [specification #40](https://github.com/Moke110/arc/issues/40)
under [ADR 0042](../../docs/adr/0042-tenv-model-benchmark-redesign.md).
The requested spelling `03_modes_tenv` is intentional. Historical benchmark
artifacts are sealed under `stage2/_archive/03_models_cmbench/`; their evidence
and the `_frozen/cm-data-v1/` bridge are retained.

The input release projects tenv raw counts to the ordered CM gene axis before
CP10K/log1p normalization. Same-plate DMSO means are allowed baseline inputs.
Training and validation targets are Δ from those controls; test-treated Δ and
directional labels live only under `fold_N/evaluator/`. Predictor manifests
contain chemical, dose, case and approved control paths only. Each group's
well predictions and truth are averaged before equally weighted group scoring.

All entrants share chirality-aware radius-2 binary Morgan-2048 vectors and three
parent-connectivity-block folds (split seed 20261105), crossed with training
seeds 17, 42 and 2026. Model admission requires separate runtime evidence on all
18,164 genes; the discovery ledger does not establish successful admission.

Build and verify inputs from the repository root using the `arc` interpreter:

```bash
python stage2/03_modes_tenv/execution/script/inputs.py build
python stage2/03_modes_tenv/execution/script/inputs.py verify
python -m unittest discover -s stage2/03_modes_tenv/execution/tests -p test_inputs.py -v
```

`output/input/manifest.json` records exact membership counts, source/config hashes,
RDKit version, structure normalization and all release file hashes. Arrays use
float32 `.npy` files for bounded memory mapping. Source corpora and concluded
pretraining/transfer results remain owned by their existing steps.

Frozen input evidence (#42): **460 drugs, 4,206 treated wells, 1,380 drug–dose
groups and 18,164 genes**. Controls separate into 192 DMSO, 80 media and 32
blank rows. Test drug counts are 154/153/153; each validation set has 46 drugs;
training counts are 260/261/261. Two complete builds have identical manifests
and array hashes. Three synthetic CLI tests pass, including normalization
order and malformed-control rejection. Build peaks and commands are in
`execution/input_record/verification.json`; arrays remain local-only.

Resource-owned execution uses the isolated artifact adapters through one scheduler:

```bash
/home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/schedule.py run <queue.json>
```

A queue records `input_root`, `output_root` and unique model/fold/seed `jobs`.
Set `max_active_jobs: 1` for long training runs; the default ceiling is two,
subject to real host/GPU telemetry, measured resource profiles and headroom.
Unknown profiles start alone. Each job may set `inference_timeout_seconds`
(default 1800); inference has its own scheduler watchdog and does not consume
its cumulative four-hour training allowance. Terminal recovery outcomes are
published to both scheduler state and `run.json`, so dead or interrupted workers
remain visible to the prediction barrier. Completed outputs still require
identity, prediction and checkpoint verification before reuse.


The runtime entrant registry is frozen in
[`execution/formal/registry.json`](execution/formal/registry.json): PRnet, CPA,
chemCPA, CMonge and PerturbDiff passed real full-axis train/save/reload/predict
and resource admission. Eleven candidates remain excluded with evidence and
four remain unresolved with named gaps in that registry. This admission proves
execution under the selected contract; it does not establish convergence or
test performance. Native setup evidence remains distinct from runtime admission.

The formal matrix contains **63 expected outcomes**: five predictors and the
`no_change`/`delta_mean` references, each with three folds and three seeds. The
reference seed copies are deterministic. Commands, resource profiles, input
hashes and expected identities are pinned under `execution/formal/`; the queue
uses serial long training and each run's four-hour cumulative training budget.
The native measured recipes imply roughly two days of serial work for CMonge
and PerturbDiff alone. These projections are estimates; actual terminal statuses
and telemetry determine completion.

Historical v1 dispatch command (joint v2 migration is now applied; use the
reviewed joint runtime for the next dispatch, after its separate preflight):

```bash
/home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/schedule.py run \
  stage2/03_modes_tenv/execution/formal/queue.json
```

Predictions and checkpoints stay local at `output/formal/`. Canonical final test
scoring waits for all 63 completed or terminal failure outcomes; no pilot is a
formal run. The user-authorized 2026-11-06 amendment in ADR 0042 permits separate
partial scoring of validated completed predictions while remaining fits continue.
After that barrier, use `benchmark.py seal` with the frozen commands and models,
then `evaluate --workers 2` and `replay`. All three figures export both PNG and
text-preserving SVG, with metric-defined counts distinguished from completed
runs.


For an unsupplemented evaluation, publish the verified CSV tables and three
PNG/SVG figures to the step output directory:

```bash
/home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/deliver.py \
  --output stage2/03_modes_tenv/output/formal \
  --publish stage2/03_modes_tenv/output
```

This artifact-only delivery replays the original checked evaluation, then adds
completion, incomplete-entry and deterministic-reference labels to the best-run
figure. It preserves all scored values, selected runs, source arrays and frozen
run identities. `publication.json` records source/evaluation/rendering hashes
and the exact published table/figure hashes. Run delivery again after an original
`benchmark.py replay` to restore the additional presentation labels. Canonical final publication follows the formal evaluation barrier. Interim
scores, provenance and snapshots use a separate partial-evaluation namespace.
User-facing interim figures are `output/metrics.png`, `output/pr_panels.png` and
`output/best_run.png`, with SVG companions when produced. The canonical final
figures overwrite those same display paths after all 63 outcomes settle; the
partial snapshot evidence remains separate.

Native learning curves are an additional artifact-only deliverable. Rebuild the
interim figure without loading models, datasets, prediction arrays or test truth:

```bash
/home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/learning_curves.py
```

For unsupplemented final delivery and replay, run that command with `--final`
after all 63 outcomes are terminal. The supplement-aware replay below refreshes
the final learning curves alongside the other figures. The renderer checks the
terminal barrier and writes `output/learning_curves.png` and its text-preserving
SVG companion at the same display paths. `output/learning_curves.json` records source checkpoint
hashes, actual plotted coordinates, available/missing histories, run statuses and
the renderer hash. These files are presentation artifacts, not formal scores.
The renderer is standalone and outside the frozen common-source and adapter
identity allowlists; training code, recipes and existing replay code stay frozen.

Each fitted-model panel expects nine blue training and nine green validation
curves (three folds × seeds 17/42/2026); available counts are printed in-panel.
PRnet uses native Gaussian NLL; CPA uses its saved reconstruction NLL component;
chemCPA uses its saved combined autoencoder objective; CMonge uses the native
Sinkhorn-plus-Monge objective; PerturbDiff uses its saved weighted diffusion loss.
The axes preserve native epoch/iteration/step coordinates and objective units.
Native reductions differ, so values are not comparable across models. Formal CMonge v2 averages all 138 validation conditions with fixed draws;
legacy v1 samples one held-out condition per event. Formal PerturbDiff v2 pairs
its last training batch with all 138 group-weighted EMA validation groups; legacy
v1 covers the first 40 groups. Both models have separate legacy panels excluded
from their nine formal slots. CMonge durable `history.json` can show events before
final checkpoint/inference; source hashes and native coordinates remain explicit. Histories can remain unavailable until an
adapter completes training and inference. Missing/unfinished histories are
omitted without smoothing or reconstruction. `no_change` and `delta_mean` have
explicit “not applicable” panels because deterministic references do not fit.


Historical continuation **resumed at the user's request on 2026-11-06** from the verified
34/63-run save point, using the unchanged serial queue. Canonical final test scoring
remains withheld until all 63 outcomes settle; partial evaluation is permitted
under the ADR 0042 amendment. The historical pause occurred after 34/63 formal
runs completed naturally (both references 9 each, PRnet 9, CPA 7), with zero failures
and 29 not dispatched. No formal test scoring has occurred. The verified pause
manifest, remaining cases and exact same-queue continuation commands are in
[`execution/formal/PAUSE.md`](execution/formal/PAUSE.md). Full binaries remain
local; the pause record describes the stopped state before this continuation.


CMonge and PerturbDiff are each assigned nine formal v2 slots under the
[ADR 0042 joint amendment](../../docs/adr/0042-tenv-model-benchmark-redesign.md#amendment-2026-11-06-joint-cmonge-and-perturbdiff-v2-protocols).
The joint 18-identity migration is **APPLIED** after independent review and
protected atomic exchange. Fresh postproof verified all 63 identities (18 changed,
45 unchanged), the unchanged 45 queue records, all 22 archive operations and both
formal inventories. The original 14,400-second cumulative caps survive in the
tracked `execution/formal/budget-ledger.json`; the tracked migration plan records
the reviewed exchange. Retained `execution/.formal-v2-generation/` recovery bytes
and native legacy binaries remain local. Four completed legacy outcomes and the
interrupted PerturbDiff outcome supply no formal v2 completion credit.

At the 2026-11-06 23:29:54 actual-host check, the owned scheduler service was
**ACTIVE** (MainPID 1487683, `Restart=no`) in runtime startup validation; first-worker
dispatch/account/telemetry proof remained pending. This is service activation,
not v2 completion. The earlier automatic approval rejection preceded explicit
user confirmation; the joint APPLY passed automatic review. The journal and
reviewed inventory hashes are recorded in ADR 0042.

The initial startup subsequently failed closed before worker dispatch because
`nvidia-smi` was absent from the service PATH. Its first budget account had no
attempts and retained the historical debit unchanged. The launch owner verified
the old service inactive/MainPID 0 and lock free, then restarted the same owned
unit with `/usr/lib/wsl/lib` on PATH (invocation
`9cc05786bff54bb3a037840a21316c30`). Integrity checks were running again;
first-worker/account/telemetry proof remained pending at this documentation
snapshot. No protocol, formal identity or cumulative cap changed for that fix.

CMonge validates all 138 equal-weight conditions with fixed draws and saves
latest/best prediction-only checkpoints and event history every 50 iterations.
Minimum 500 iterations, patience 10 subsequent checks and relative delta 0.001
give earliest patience stop at 1,000 iterations. PerturbDiff validates all 138
group-weighted EMA groups and saves checkpoints every 1,000 steps; minimum
20,000, patience 10 and absolute delta 0.01 give earliest stop at 30,000 steps.
Both restart from seed under the original cumulative 14,400-second allowance,
including legacy attempts and retries; exact resume is unsupported. CMonge's
measured full-validation maximum projects ~6.14 hours, and fold 0 seed 17 has
only 424.49 seconds remaining. Preserve full validation and report honest
terminal timeouts under independently enforced wall deadlines.

The standalone learning-curve renderer supports both v2 native histories and
separate legacy panels (commit `d03f44a`); its existing curve contract above
applies. Missing histories stay explicit and legacy trajectories are never
spliced into formal v2 slots. Pre-migration legacy statuses do not satisfy the
final 63-outcome barrier. #40/#46/#47 stay open through formal delivery/review.

## CMonge v2 fold 0 seed 17 supplemental attempt

The [ADR 0042 supplemental amendment](../../docs/adr/0042-tenv-model-benchmark-redesign.md#amendment-2026-11-07--cmonge-v2-fold-0-seed-17-supplemental-attempt)
authorizes a single fresh **14,400-second training allowance** and a separate
**1,800-second inference allowance**, using the same CMonge v2 recipe from
seed initialization. The original cumulative-budget timeout and its artifacts
remain in `output/formal/`. The supplemental runner
[`cmonge_supplement.py`](execution/script/cmonge_supplement.py) uses a separate
attempt directory and durable claim.

[`supplement_overlay.py`](execution/script/supplement_overlay.py) binds the
reviewed supplemental outcome to the existing fold/seed slot and the other
62 original outcomes. An authenticated timeout or failure occupies the slot
without a score. Invalid evidence blocks final evaluation. Run supplement-aware
publication after scoring and on every replay to retain the budget exception,
measured phase times and host RSS in the final tables and figures. The
learning curves display the original timeout separately from the selected
supplement and the legacy histories.

### Final results

The original scheduler exited cleanly with 62 completions and the CMonge
timeout. The supplemental attempt then stopped by validation patience at
iteration 2,150 after 43 validation events with 138 conditions per event.
Watchdog setup and training took 4,253.63 s, including 4,228.00 s of native
training. Observed inference took 1.19 s, including 0.72 s of native inference.
Peak host RSS was 1.68 GiB and sampled peak workstation GPU usage was 2.80 GiB.
The original timeout's host RSS is unavailable. Its history has four validation
events and is shown separately in the learning curves.

The selected matrix completed **63/63** slots across five admitted predictors
and two deterministic references. The candidate registry records eleven
excluded and four unresolved candidates. The evaluator produced 28,980
drug–dose group rows and 63 run-level PR sources. Every model–metric pair has
nine defined runs except `no_change` PCC and Spearman, which are undefined in
all nine runs because its predictions are constant. Each reference has three
folds with deterministic seed copies.

PerturbDiff led seven of the eight nine-run metric means, including positive
and negative AP@100 of 0.23914 and 0.23463 and RMSE of 0.13172. The training
group mean reference led mean R² at 0.02252. CMonge's directional AP@100 means
were 0.18588 and 0.18040 under the disclosed supplemental budget. The best-run
figure is a post-hoc view. The summary describes held-out performance on the
frozen tenv release, with variability across folds and seeds.

- Metrics [PNG](output/metrics.png) and [SVG](output/metrics.svg)
- Directional PR panels [PNG](output/pr_panels.png) and [SVG](output/pr_panels.svg)
- Best-run PR [PNG](output/best_run.png) and [SVG](output/best_run.svg)
- Native learning curves [PNG](output/learning_curves.png), [SVG](output/learning_curves.svg) and [coordinates](output/learning_curves.json)
- Tables for [groups](output/per_group.csv), [runs](output/per_run.csv), [summaries](output/summary.csv), [PR bands](output/pr_bands.csv) and [best runs](output/best_runs.csv)
- [Publication hashes](output/publication.json), [slot provenance](output/slot_provenance.json) and [sealed prediction manifest](execution/final_record/predictions.json)

The original joint v2 scheduler sampled host pressure and the read-only watcher
monitored it. Automatic in-flight pressure stops were unavailable in that
scheduler, and its resource admission check preceded input verification.
The original watchdog charged full worker wall, including inference, against
the remaining training allowance while also applying an inference sublimit.
Its terminal `training_seconds` is the cumulative debit. Per-attempt GPU peaks
are absent from those terminal records. The supplement has separate phase
deadlines and resource evidence. The preserved original account, run records,
legacy archive and corpus release are read-only inputs to final delivery.
Scoring and replay write into the private overlay and presentation paths.

### Replay the published figures

Use the preserved workstation artifacts, with the scheduler idle and at least
12 GiB available RAM. Run from the repo root. The pinned publication command
uses the saved evaluator results and refreshes all four figure pairs and their
provenance. The initial run used `prepare`, `verify` and `evaluate` in order,
with independent review of the emitted selection and evaluation digests.

```bash
step=/home/chang/arc/stage2/03_modes_tenv
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- \
  env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
  /usr/bin/time -v /home/chang/miniconda3/envs/arc/bin/python \
  "$step/execution/script/supplement_overlay.py" publish \
  --formal "$step/execution/formal" \
  --reviewed-stage-sha256 73fb21036a2d04a12c5d8d9ce712ed544b1c38eb70fbfff884dc6753c8b60028 \
  --policy "$step/execution/final_record/policy.json" \
  --reviewed-policy-sha256 3b577b7a441b6e0f730c808b79b98ef57bd79c6a13a3b2dfa9e9ff29550923db \
  --reviewed-terminal-sha256 591872e3aca0a6fba30b0b308381a8569b62c24b1a6ade07a7289c69f0eefb67 \
  --reviewed-selection-sha256 e2bc2ea2414e52a905c405d37f7ff58acab48e6c56211fd0f4c9d1269bbbae56 \
  --reviewed-evaluation-sha256 839045cc974386d1c6d21c61753ede02e92a0b208845816f090a94e903f6a918
```

Preparation, two-worker evaluation and publication peaked at 210,896,
422,864 and 259,312 KiB RSS, respectively. The [final evidence index](execution/final_record/evidence.sha256)
binds the compact audit records and resource logs. The live prediction arrays,
model weights and PR source array remain local.
