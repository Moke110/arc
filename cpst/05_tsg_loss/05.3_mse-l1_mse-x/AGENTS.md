# Step 05.3 — TranSiGen reconstruction MSE sweep

Read the [study report](README.md), [Step 05 routing](../AGENTS.md) and
[ADR 0091](../../../docs/adr/0091-transigen-reconstruction-mse-sweep.md).
For runner or scoring changes, use the
[shared TranSiGen APIs](../../util/README.md#transigen-training-and-evaluation).

## Inputs and protocol

Use the frozen [Denoise-DEG 5k release](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json),
[tenv cold-drug folds 0–2](../../03_data_train/output/frozen-folds/tenv/manifest.json),
seed 42, frozen KPGT vectors and
[TranSiGen vendor source](../../../pert_model/TranSiGen/manifest.sha256).
Use measured same-plate DMSO as the model input and the reference for both
Delta vectors. Score with the frozen HVG-DEG 5k training bandwidths.
[Step 05.2](../05.2_mse_l1/AGENTS.md) documents the inherited protocol.

Use the [registered weight grid](execution/sweep_protocol.py) for the shared
control and treatment reconstruction MSE weight `alpha`. Keep Delta MSE at
unit weight, Delta L1 at 1.75 and the summed native KL terms at 0.1. At
`alpha = 1`, preserve Step 05.2 objective arithmetic exactly. Select the
earliest minimum per-condition registered validation objective. Pair held-out
scores by fold and seed against the fresh unit-weight fit.

## Prepare the runtime

Activate `arc` at the materialized checkout root and define the
[Step 04 capped runtime](../../04_tsg_axis/AGENTS.md#prepare-the-capped-runtime).
Follow the [execution environment rules](../../../docs/agents/execution-environment.md).
Keep candidates under `tmp/`, run one CUDA fit at a time and record WSL RAM
and GPU usage every ten seconds.

Pass `--training`, `--representations` and `--source` to the runner, matrix
or formal report when the pinned inputs are materialized elsewhere. The evaluator
and formal report accept `--training` and `--bandwidth-training`. Use
`--previous-step` on the formal report to locate the accepted Step 05.2 root.

## Run the bounded comparison

Use [run_reconstruction.py](execution/run_reconstruction.py) for two-epoch
fits at weights 0 and 1 on fold 0. Register both, then run `train`, `predict`
and `replay` in separate processes for each coordinate. Use
[evaluate_reconstruction.py](execution/evaluate_reconstruction.py) to score
and verify the fits. Compare two renders from the smoke report CSVs.
[Smoke receipt](log/335-smoke.json)

```bash
(
  set -e
  step="$PWD/cpst/05_tsg_loss/05.3_mse-l1_mse-x"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-reconstruction-smoke.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib" "$work/csv-only"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
  runs=()
  for alpha in 0 1; do
    run="$work/alpha-$alpha"
    capped_step04 python "$step/execution/run_reconstruction.py" register \
      --alpha "$alpha" --fold 0 --smoke --output "$run"
    runs+=(--run "$run")
  done
  for alpha in 0 1; do
    for phase in train predict replay; do
      capped_step04 python "$step/execution/run_reconstruction.py" "$phase" \
        --output "$work/alpha-$alpha"
    done
  done
  for phase in score verify; do
    capped_step04 python "$step/execution/evaluate_reconstruction.py" "$phase" \
      --smoke "${runs[@]}" --output "$work/report" --scratch "$work/scratch"
  done
  cp "$work/report/"{delta-metrics,epoch-objectives}.csv "$work/csv-only/"
  for attempt in first second; do
    for renderer in render_metrics render_learning_curves; do
      python "$step/execution/$renderer.py" --smoke \
        --data "$work/csv-only" --output "$work/$attempt"
    done
  done
  for name in delta-metrics.png delta-metrics-improvement.png learning-curves.png; do
    cmp "$work/first/$name" "$work/second/$name"
  done
)
```

## Reproduce the original formal sweep

[matrix_reconstruction.py](execution/matrix_reconstruction.py) registers all
27 coordinates before fitting, checks sealed phases when resuming at the same
matrix and audit paths, and records peak RSS and ten-second WSL RAM and GPU
samples. Allow at least 4 GiB of free space. Set `TSG_RECONSTRUCTION_WORK_ROOT`
to an existing directory or backed symlink under checkout `tmp/` to use a
different volume. Keep scratch, report and publication on the same filesystem.

[report_reconstruction.py](execution/report_reconstruction.py) checks exact
independent replay, nine held-out metrics per coordinate with eligibility
counts, frozen bandwidths and byte-identical CSV-only rendering. It matches
the three fresh unit-weight fits to [sealed Step 05.2](../05.2_mse_l1/log/333-formal.json)
at Delta L1 weight 1.75 on traces, selected and stopping epochs, validation
and test predictions, Delta arrays, scores and eligibility. The build also
includes the coordinator resource capture if the same run's
`coordinator-resource-samples.jsonl` is present in the matrix audit directory.
[Formal receipt](log/336-formal.json)

```bash
(
  set -e
  step="$PWD/cpst/05_tsg_loss/05.3_mse-l1_mse-x"
  mkdir -p "$PWD/tmp"
  work_root="${TSG_RECONSTRUCTION_WORK_ROOT:-$PWD/tmp}"
  work=$(mktemp -d "$work_root/tsg-reconstruction-formal.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
  revision=$(git rev-parse HEAD)
  capped_step04 python -u "$step/execution/matrix_reconstruction.py" \
    --output "$work/matrix" --logs "$work/matrix-audit"
  report_args=(
    --matrix "$work/matrix" --matrix-audit "$work/matrix-audit"
    --report "$work/report" --publication "$work/publication"
    --activity "$work/report-audit" --scratch "$work/scratch"
  )
  capped_step04 python "$step/execution/report_reconstruction.py" build \
    "${report_args[@]}" --code-revision "$revision"
  capped_step04 python "$step/execution/report_reconstruction.py" verify \
    "${report_args[@]}"
)
```

## Replay a sealed fit

Use the registered environment and original `stage2/` paths at revision
`fa7e9a0e580d06a63210e93bff82c655d37b5ae8`. Replay a scratch copy and keep
installed evidence read-only. This example uses the unit-weight fit on fold 0.
[Registration verification](../../util/training/verification.md)

```bash
(
  set -e
  stage="$PWD/cpst"
  study="05_tsg_loss/05.3_mse-l1_mse-x"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-reconstruction-replay.XXXXXX")
  git worktree add --detach "$work/code" fa7e9a0e580d06a63210e93bff82c655d37b5ae8
  trap 'git worktree remove "$work/code"' EXIT
  mkdir "$work/scratch" "$work/run"
  cp --reflink=auto --dereference -- \
    "$stage/$study/output/formal/matrix/runs/alpha-1/fold-0/"* "$work/run/"
  export TMPDIR="$work/scratch" PYTHONDONTWRITEBYTECODE=1
  capped_step04 python "$work/code/stage2/$study/execution/run_reconstruction.py" replay \
    --output "$work/run" \
    --training "$stage/03_data_train/03.7_tenv_denoise-5k/output" \
    --representations "$stage/util/output/representations" \
    --source "$PWD/pert_model/TranSiGen/source"
)
```

## Verify the current figures

Render twice from the released plot CSVs and compare both renders to the
accepted PNGs. [CSV-only renderer](figure/script/render_comparison.py) ·
[Render receipt](log/figure-render.json)

```bash
(
  set -e
  figure="$PWD/cpst/05_tsg_loss/05.3_mse-l1_mse-x/figure"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-reconstruction-figures.XXXXXX")
  export TMPDIR="$work"
  for attempt in first second; do
    MPLCONFIGDIR="$work/mpl-$attempt" python "$figure/script/render_comparison.py" \
      --data "$figure/data" --output "$work/$attempt"
    for name in learning-curve.png metrics.png metrics-paired.png; do
      cmp "$figure/$name" "$work/$attempt/$name"
    done
  done
)
```

## Rebuild the current figures

Use the [shared comparison builder](../execution/build_comparison.py) to
rebuild the plot CSVs and PNGs from the accepted fits. Pass
`--evidence-root` for a different checkout with the accepted inputs.

```bash
(
  set -e
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-reconstruction-figure-build.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
  capped_step04 python cpst/05_tsg_loss/execution/build_comparison.py \
    --step 05.3_mse-l1_mse-x --output "$work/candidate"
)
```
