# Step 06.4 TranSiGen measured encoder targets

Read the [Round 4 design](../../../docs/adr/0098-transigen-shape-round4.md)
and [inherited Round 3 protocol](../../../docs/adr/0097-transigen-shape-round3.md)
before changing an arm, objective or evaluation. Use the shared
[TranSiGen training and evaluation API](../../util/README.md#transigen-training-and-evaluation)
through the numbered entry points. The [study report](README.md) owns the
accepted comparisons and evidence.

## Register and fit

1. Use the [tenv Denoise-DEG 5k release](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json),
   frozen tenv KPGT vectors from the [representation release](../../util/README.md),
   and [author source](../../../pert_model/TranSiGen/manifest.sha256).
   Keep the cold-drug folds, measured same-plate DMSO and training-fold dose
   transform fixed. Use the [tenv HVG-DEG 5k release](../../03_data_train/03.1_tenv_hvg-5k/output/manifest.json)
   for fold bandwidth calibration.
2. Follow the [bounded CUDA gate](#bounded-cuda-gate) before the
   formal matrix. Use one capped CUDA coordinate at a time, with WSL RAM,
   Windows host RAM and GPU use sampled every ten seconds. Follow the
   [execution environment rules](../../../docs/agents/execution-environment.md).
3. Use [the matrix runner](execution/matrix_measured.py) to register the
   complete fit roster before training. It calls
   [the coordinate runner](execution/run_measured.py) for `register`, `train`,
   `predict` and `replay` in separate processes. Keep training seeds separate
   from the frozen validation and test sampling seeds. In direct Δ arms,
   form treated predictions by adding decoded Δ to measured control.

## Prepare a bounded or formal run

Run from the materialized repository root with `arc` active. Define these
functions in the shell where you will run the gate or formal build. Keep
the accepted Step 06.1–06.3 releases read-only for lineage verification.

```bash
prepare_round4() {
  STEP="$PWD/cpst/06_tsg_module/06.4_round4"
  mkdir -p "$PWD/tmp"
  WORK=$(mktemp -d "$PWD/tmp/tsg-measured-$1.XXXXXX")
  mkdir -p "$WORK/scratch" "$WORK/matplotlib"
  export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib" PYTHONDONTWRITEBYTECODE=1
}

capped_round4() {
  LC_ALL=C free -g | awk '$1 == "Mem:" {available=$7} END {exit available<12}' || return
  systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
    env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 "$@"
}

stop_sampler() {
  local status=0
  kill -TERM "$SAMPLER_PID" 2>/dev/null || true
  wait "$SAMPLER_PID" || status=$?
  trap - EXIT
  return "$status"
}

start_sampler() {
  python "$STEP/../06.1_round1/execution/sample_resources.py" \
    --output "$WORK/coordinator-resources.jsonl" &
  SAMPLER_PID=$!
  trap 'stop_sampler || true' EXIT
  until [ -s "$WORK/coordinator-resources.jsonl" ]; do
    kill -0 "$SAMPLER_PID" || return 1
    sleep 1
  done
}
```

## Bounded CUDA gate

Require the smoke matrix, objective audit, independent score verification,
prior lineage check, paired table and figure checks, model tests and
resource audit to pass before a formal build. Use the
[gate at `f387f2c3`](https://github.com/Moke110/carps/blob/f387f2c370eb87cd27c02502aa98f165851491a7/stage2/06_tsg_shape/06.4_round4/README.md#bounded-cuda-gate)
for the renderer and command versions that produced the
[accepted gate hashes](log/358-gate.json).

```bash
(
  set -e
  prepare_round4 smoke
  start_sampler
  capped_round4 python "$STEP/execution/matrix_measured.py" --smoke \
    --output "$WORK/matrix" --logs "$WORK/matrix-audit"
  sleep 10
  capped_round4 python "$STEP/execution/audit_measured.py" objectives --smoke --matrix "$WORK/matrix"
  for phase in score verify; do
    capped_round4 python "$STEP/execution/evaluate_measured.py" "$phase" --smoke \
      --matrix "$WORK/matrix" --output "$WORK/report" --scratch "$WORK/scratch"
  done
  capped_round4 python "$STEP/execution/prior_measured.py"
  python "$STEP/execution/audit_measured_scores.py" --smoke --report "$WORK/report"
  python "$STEP/execution/audit_measured_figures.py" build --smoke \
    --data "$WORK/report" --output "$WORK/figures" --receipt "$WORK/figures.json" \
    --scratch "$WORK/scratch"
  capped_round4 env TSG_MEASURED_TEST_DEVICE=cuda python -m unittest discover \
    -s "$STEP/../../util/training/tests/06.4_round4" -p 'test_measured_model.py'
  stop_sampler
  python "$STEP/execution/audit_measured.py" resources --smoke \
    --matrix-audit "$WORK/matrix-audit" --coordinator "$WORK/coordinator-resources.jsonl" \
    --coverage-output "$WORK/coordinator-coverage.jsonl"
)
```

## Reproduce the formal study

For new registrations, verify the 18 retained fits with the
[lineage verifier](execution/prior_measured.py) before fitting. Use the
[publisher](execution/publish_measured.py) to run the
[evaluator](execution/evaluate_measured.py) in separate score and verify
processes. Require 27 exact replays, independent agreement on all 405
held-out metric values and a paired contrast audit within fold and seed.

```bash
(
  set -e
  prepare_round4 formal
  capped_round4 python "$STEP/execution/prior_measured.py" > "$WORK/source-lineage.json"
  start_sampler
  capped_round4 python "$STEP/execution/matrix_measured.py" \
    --output "$WORK/matrix" --logs "$WORK/matrix-audit"
  sleep 10
  stop_sampler
  capped_round4 python "$STEP/execution/publish_measured.py" build \
    --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
    --report "$WORK/report" --publication "$WORK/publication" \
    --activity "$WORK/publication-audit" --scratch "$WORK/scratch" \
    --coordinator-samples "$WORK/coordinator-resources.jsonl" --code-revision "$(git rev-parse HEAD)"
)
```

To verify [the installed #359 release](log/359-formal.json), use `verify`
with the [original publisher and renderer at `09da1ead`](https://github.com/Moke110/carps/blob/09da1ead4ee31d5c8d84ac60d252a4a2f3cc209f/stage2/06_tsg_shape/06.4_round4/README.md#reproduce-the-formal-study)
from that materialized checkout with the original `stage2/` layout.
The historical verifier reads the installed results, independently
rescores, audits paired tables and renders the original figures twice.

```bash
(
  set -e
  prepare_round4 verify
  STEP="$PWD/stage2/06_tsg_shape/06.4_round4"
  FORMAL="$STEP/output/formal"
  capped_round4 python "$STEP/execution/publish_measured.py" verify \
    --matrix "$FORMAL/matrix" --matrix-audit "$FORMAL/matrix-audit" \
    --report "$FORMAL/report" --publication "$STEP" \
    --activity "$WORK/publication-audit" --scratch "$WORK/scratch"
)
```

## Replay a sealed Step 06.4 fit

The [formal receipt](log/359-formal.json) pins the registered fit code at
`f387f2c370eb87cd27c02502aa98f165851491a7`. The registration binds the
original `stage2/` paths and code digests. From the materialized current
checkout with `arc` active, copy the shared-decoder fit for fold 0,
seed 17 and give the original runner the current locations of its inputs.

```bash
(
  set -e
  BASE="$PWD"
  mkdir -p "$BASE/tmp"
  WORK=$(mktemp -d "$BASE/tmp/tsg-measured-pinned.XXXXXX")
  STAGE="$BASE/cpst"
  STEP="$STAGE/06_tsg_module/06.4_round4"
  SEALED_RUN="$STEP/output/formal/matrix/runs/XΔ-XΔ/fold-0/seed-17"
  git worktree add --detach "$WORK/code" f387f2c370eb87cd27c02502aa98f165851491a7
  cleanup_replay() {
    if [ -n "${REPLAY_SAMPLER_PID:-}" ]; then
      kill -TERM "$REPLAY_SAMPLER_PID" 2>/dev/null || true
      wait "$REPLAY_SAMPLER_PID" || true
    fi
    git worktree remove "$WORK/code"
  }
  REPLAY_SAMPLER_PID=
  trap cleanup_replay EXIT
  mkdir -p "$WORK/scratch" "$WORK/run"
  cp --reflink=auto --dereference -- "$SEALED_RUN"/* "$WORK/run/"
  export TMPDIR="$WORK/scratch" PYTHONDONTWRITEBYTECODE=1

  python "$STEP/../06.1_round1/execution/sample_resources.py" \
    --output "$WORK/coordinator-resources.jsonl" &
  REPLAY_SAMPLER_PID=$!
  until [ -s "$WORK/coordinator-resources.jsonl" ]; do
    kill -0 "$REPLAY_SAMPLER_PID" || exit 1
    sleep 1
  done
  LC_ALL=C free -g | awk '$1 == "Mem:" {available=$7} END {exit available<12}'
  systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
    env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
      python "$WORK/code/stage2/06_tsg_shape/06.4_round4/execution/run_measured.py" replay \
      --output "$WORK/run" \
      --training "$STAGE/03_data_train/03.7_tenv_denoise-5k/output" \
      --representations "$STAGE/util/output/representations" \
      --source "$BASE/pert_model/TranSiGen/source"
  kill -0 "$REPLAY_SAMPLER_PID"
)
```

## Replay the published figures

Use the `arc` environment and the package versions in the
[figure receipt](log/figure-render.json). The current renderer uses the
[shared plotting API](../../util/README.md#plotting) and the three plotted
CSVs. Compare all three PNGs byte for byte with the published images.

```bash
(
  set -e
  STEP="$PWD/cpst/06_tsg_module/06.4_round4"
  mkdir -p "$PWD/tmp"
  WORK=$(mktemp -d "$PWD/tmp/tsg-measured-figures.XXXXXX")
  mkdir -p "$WORK/scratch" "$WORK/matplotlib"
  export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib" PYTHONDONTWRITEBYTECODE=1
  python "$STEP/figure/script/render_comparison.py" \
    --data "$STEP/figure/data" --output "$WORK/figures"
  for figure in learning-curve.png metrics.png metrics-paired.png; do
    cmp "$STEP/figure/$figure" "$WORK/figures/$figure"
  done
)
```
