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

Dispatch from the repository root after offline gates pass:

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


After evaluation, publish the verified CSV tables and three PNG/SVG figures to
the step output directory (user-selected display location):

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


Execution **resumed at the user's request on 2026-11-06** from the verified
34/63-run save point, using the unchanged serial queue. Canonical final test scoring
remains withheld until all 63 outcomes settle; partial evaluation is permitted
under the ADR 0042 amendment. The historical pause occurred after 34/63 formal
runs completed naturally (both references 9 each, PRnet 9, CPA 7), with zero failures
and 29 not dispatched. No formal test scoring has occurred. The verified pause
manifest, remaining cases and exact same-queue continuation commands are in
[`execution/formal/PAUSE.md`](execution/formal/PAUSE.md). Full binaries remain
local; the pause record describes the stopped state before this continuation.
