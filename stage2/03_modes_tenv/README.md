# 03_modes_tenv — tenv shared-Morgan comparison

## Retrieval evaluation

The active retrieval statistics are **AUPRC**, **AUPRC_positive**, and
**AUPRC_negative**, with corresponding **PR curve** panels. All use full-CM
DESeq2 significance and direction matching, under [ADR 0059](../../docs/adr/0059-stage2-significance-retrieval.md).
The [shared evaluator](../evaluation/README.md) defines the denominator,
undefined groups, pooling, and bounded resource policy.

- [PR curve](output/retrieval/PR_curve.png)
- [AUPRC and continuous metrics](output/retrieval/AUPRC.png)
- [Run statistics](output/retrieval/per_run.csv), [dose statistics](output/retrieval/per_dose.csv), and [summary](output/retrieval/summary.csv)
- [Replay manifest](output/retrieval/manifest.json)

From the repository root, after the relevant significance release is ready:

```bash
/home/chang/miniconda3/envs/arc/bin/python stage2/evaluation/publish.py --step 3
/home/chang/miniconda3/envs/arc/bin/python stage2/evaluation/evaluate.py replay --step 3
```

Previous retrieval scores, metric figures, and historical copies were removed;
Git retains their version history. Existing checkpoints and their original
validation-loss selections are unchanged. Other sections below describe input
and training provenance, rather than the retired retrieval publication.

This step implements [specification #40](https://github.com/Moke110/arc/issues/40)
under [ADR 0042](../../docs/adr/0042-tenv-model-benchmark-redesign.md).
The requested spelling `03_modes_tenv` is intentional. The prior benchmark's
tracked records are available in Git at `f95fe1e`. Its archival manifest remains
in `execution/input_record/archive.json`.

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
float32 `.npy` files for bounded memory mapping. Source corpora remain owned by `02_dataset_corpus/`.

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

Predictions and checkpoints stay local at `output/formal/`. The original
63-slot terminal barrier and authenticated completed-run verification remain
available. Current tenv scoring, figures and replay use the shared commands
above. Earlier partial/final score publication and best-run figures are retired.

Native learning curves are an additional artifact-only deliverable. Rebuild the
interim figure without loading models, datasets, prediction arrays or test truth:

```bash
/home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/learning_curves.py
```

For unsupplemented final delivery and replay, run that command with `--final`
after all 63 outcomes are terminal. The standalone renderer refreshes the final learning curves. The renderer checks the
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
omitted without smoothing or reconstruction. The learning figure includes only
fitted models and their available legacy histories.


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
learning figure plots the selected supplemental history alongside the other
formal histories and retains separate legacy panels. The original timeout
history remains in the JSON evidence and is omitted from the figure.

### Final results

The original scheduler exited cleanly with 62 completions and the CMonge
timeout. The supplemental attempt then stopped by validation patience at
iteration 2,150 after 43 validation events with 138 conditions per event.
Watchdog setup and training took 4,253.63 s, including 4,228.00 s of native
training. Observed inference took 1.19 s, including 0.72 s of native inference.
Peak host RSS was 1.68 GiB and sampled peak workstation GPU usage was 2.80 GiB.
The original timeout's host RSS is unavailable. Its four validation events
remain in the learning-curve JSON evidence.

The selected matrix completed **63/63** slots across five admitted predictors
and two deterministic references. The candidate registry records eleven
excluded and four unresolved candidates. The shared evaluator now rescores these frozen predictions under ADR 0059.
Each reference has three folds with deterministic seed copies.

The original joint v2 scheduler sampled host pressure and the read-only watcher
monitored it. Automatic in-flight pressure stops were unavailable in that
scheduler, and its resource admission check preceded input verification.
The original watchdog charged full worker wall, including inference, against
the remaining training allowance while also applying an inference sublimit.
Its terminal `training_seconds` is the cumulative debit. Per-attempt GPU peaks
are absent from those terminal records. The supplement has separate phase
deadlines and resource evidence. The preserved original account, run records
and corpus release are read-only inputs to final delivery.
Scoring and replay write into the private overlay and presentation paths.
The reviewed supplemental runner derives the repository root from the fixed
step depth. Replay uses the reviewed source bytes and workstation layout
recorded in its policy.

### Replay the published figures

Use the shared tenv publication and replay commands at the top of this file.
Learning curves retain their standalone history-based workflow. Historical
score publications and PR source arrays are retired; trained weights, sealed
predictions and resource records are retained.
