# TranSiGen additive Delta L1 sweep

Read the [study report](README.md), [CarPS-Transcriptome routing](../../AGENTS.md),
and [ADR 0090](../../../docs/adr/0090-transigen-mse-l1-sweep.md).
For runner or scoring changes, use the [shared TranSiGen APIs](../../util/README.md#transigen-training-and-evaluation).

Use the frozen [Denoise-DEG 5k training release](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json),
KPGT vectors and [TranSiGen source](../../../pert_model/TranSiGen/manifest.sha256)
from [Step 05.1](../05.1_loss_metrics/README.md). Keep cold-drug folds 0–2 and
seed 42. Use measured same-plate DMSO as the model input and the reference for
both predicted and observed Delta. Score with the Step 05.1 metric definitions
and frozen HVG-DEG 5k bandwidths.

Sweep additive Delta L1 weight `alpha` from 0 through 2 in 0.25 increments
with unit native MSE weights and KL weight 0.1. At zero, use the native loss
arithmetic and still log the raw L1 component. Select the earliest minimum
of the registered validation objective. Pair held-out scores by fold and seed
against the fresh zero-weight fits.

## Prepare the runtime

Activate `arc` in the materialized repository root and define the
[capped runtime](../../04_tsg_axis/AGENTS.md#prepare-the-capped-runtime).
Keep candidates under `tmp/`, run one CUDA fit at a time, and record WSL RAM and
GPU usage every ten seconds. Follow the [execution environment rules](../../../docs/agents/execution-environment.md).

For relocated frozen inputs, pass `--training`, `--representations` and
`--source` to the runner, matrix and formal report as needed. Pass the frozen
HVG-DEG 5k release as `--bandwidth-training` to the evaluator or formal report,
and the Step 05.1 study root as `--previous-step` to the formal report.

## Run a bounded comparison

Use [run_sweep.py](execution/run_sweep.py) for two fold-0 fits at weights 0 and
1. Register both, then run `train`, `predict` and `replay` in separate processes
with the same output path for each coordinate. Use [evaluate_sweep.py](execution/evaluate_sweep.py)
to score and verify the replayed fits. Compare two CSV-only renders of the smoke
report. [Two-fit smoke receipt](log/332-smoke.json)

```bash
(
  set -e
  step="$PWD/cpst/05_tsg_loss/05.2_mse_l1"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-mse-l1-smoke.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
  runs=()
  for alpha in 0 1; do
    run="$work/alpha-$alpha"
    capped_step04 python "$step/execution/run_sweep.py" register \
      --alpha "$alpha" --fold 0 --smoke --output "$run"
    runs+=(--run "$run")
  done
  for alpha in 0 1; do
    for phase in train predict replay; do
      capped_step04 python "$step/execution/run_sweep.py" "$phase" \
        --output "$work/alpha-$alpha"
    done
  done
  for phase in score verify; do
    capped_step04 python "$step/execution/evaluate_sweep.py" "$phase" --smoke \
      "${runs[@]}" --output "$work/report" --scratch "$work/scratch"
  done
  for attempt in first second; do
    for renderer in render_metrics render_learning_curves; do
      python "$step/execution/$renderer.py" --smoke \
        --data "$work/report" --output "$work/$attempt"
    done
  done
  for name in delta-metrics.png delta-metrics-improvement.png learning-curves.png; do
    cmp "$work/first/$name" "$work/second/$name"
  done
)
```

## Reproduce the original formal sweep

[matrix_sweep.py](execution/matrix_sweep.py) checks ordered fold membership,
registers all 27 coordinates before fitting, and records phase peak RSS and
ten-second WSL RAM, GPU and disk samples. Resume with the same matrix and audit
directories. [report_sweep.py](execution/report_sweep.py) verifies the fit
receipts, all 243 fold scores, exact parity of the three zero-weight controls
against Step 05.1 native MSE at weight 1, and byte-identical CSV-only renders.
Keep scratch and the candidate publication on the same filesystem. To recheck
an existing candidate with the same code, rerun the `verify` command.
[Formal receipt](log/333-formal.json)

```bash
(
  set -e
  step="$PWD/cpst/05_tsg_loss/05.2_mse_l1"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-mse-l1-formal.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
  revision=$(git rev-parse HEAD)
  capped_step04 python -u "$step/execution/matrix_sweep.py" \
    --output "$work/matrix" --logs "$work/matrix-audit"
  report_args=(
    --matrix "$work/matrix" --matrix-audit "$work/matrix-audit"
    --report "$work/report" --publication "$work/publication"
    --activity "$work/report-audit" --scratch "$work/scratch"
  )
  capped_step04 python "$step/execution/report_sweep.py" build \
    "${report_args[@]}" --code-revision "$revision"
  capped_step04 python "$step/execution/report_sweep.py" verify "${report_args[@]}"
)
```

## Replay a sealed fit

The formal receipt pins revision `4d85563ea6d973f1573cfebaac182cb301e857ae`.
Use its original `stage2/` code paths and the registered environment. Supply
the current locations of the pinned inputs and replay a scratch copy of the
fit, leaving the installed release read-only. This example uses weight 1.75,
fold 0. [Registration verification](../../util/training/verification.md)

```bash
(
  set -e
  stage="$PWD/cpst"
  study="05_tsg_loss/05.2_mse_l1"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-mse-l1-replay.XXXXXX")
  git worktree add --detach "$work/code" 4d85563ea6d973f1573cfebaac182cb301e857ae
  trap 'git worktree remove "$work/code"' EXIT
  mkdir "$work/scratch" "$work/run"
  cp --reflink=auto --dereference -- \
    "$stage/$study/output/formal/matrix/runs/alpha-1.75/fold-0/"* "$work/run/"
  export TMPDIR="$work/scratch" PYTHONDONTWRITEBYTECODE=1
  capped_step04 python "$work/code/stage2/$study/execution/run_sweep.py" replay \
    --output "$work/run" \
    --training "$stage/03_data_train/03.7_tenv_denoise-5k/output" \
    --representations "$stage/util/output/representations" \
    --source "$PWD/pert_model/TranSiGen/source"
)
```

## Verify the current figures

Render twice from the released plot CSVs and compare each render with the
current PNGs. [CSV-only renderer](figure/script/render_comparison.py) ·
[Render receipt](log/figure-render.json)

```bash
(
  set -e
  figure="$PWD/cpst/05_tsg_loss/05.2_mse_l1/figure"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-mse-l1-figures.XXXXXX")
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
