# Step 05.1 — TranSiGen Delta loss metrics

Read the [Step 05 routing](../AGENTS.md), [study report](README.md), and
[accepted loss grid](../../../docs/adr/0089-transigen-delta-metric-weight.md).
For runner or scoring changes, use the [shared TranSiGen APIs](../../util/README.md#transigen-training-and-evaluation).

## Inputs and runtime

Activate `arc` in the materialized checkout root and define the
[Step 04 capped runtime](../../04_tsg_axis/AGENTS.md#prepare-the-capped-runtime).
Use the [tenv Denoise-DEG 5k release](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json),
[frozen tenv folds](../../03_data_train/output/frozen-folds/tenv/manifest.json),
frozen KPGT vectors and [pinned TranSiGen source](../../../pert_model/TranSiGen/manifest.sha256).
Keep the [Step 04 matched-plate protocol](../../../docs/adr/0086-tenv-axis-matched-plate-control.md),
architecture and fixed evaluation seeds. Use the measured same-plate control
as model input and as the common reference for both Delta vectors. Reserve
held-out treatment measurements for scoring. Select each checkpoint by its
registered validation objective. Keep formal inputs read-only and work under
the checkout's `tmp/`.

## Run a bounded comparison

Register the four fold-0 coordinates before fitting. Run the three phases in
separate processes with one CUDA fit at a time. [Smoke receipt](log/329-smoke.json)

```bash
(
  set -e
  step="$PWD/cpst/05_tsg_loss/05.1_loss_metrics"
  work=$(mktemp -d "$PWD/tmp/tsg-loss-smoke.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
  runs=()
  for loss in mse cos pcc l1; do
    run="$work/$loss"
    capped_step04 python "$step/execution/run_loss.py" register \
      --loss "$loss" --alpha 1 --fold 0 --smoke --output "$run"
    runs+=(--run "$run")
  done
  for loss in mse cos pcc l1; do
    for phase in train predict replay; do
      capped_step04 python "$step/execution/run_loss.py" "$phase" --output "$work/$loss"
    done
  done
  for phase in score verify; do
    capped_step04 python "$step/execution/evaluate_loss.py" "$phase" --smoke \
      "${runs[@]}" --output "$work/report" --scratch "$work/scratch"
  done
)
```

For relocated inputs, pass `--training`, `--representations` and `--source`
to the runner as applicable. The evaluator accepts `--training` and
`--bandwidth-training` for the corresponding releases.

## Reproduce the original formal matrix

The current numbered commands reproduce the original study design in a new
candidate. The matrix driver validates all 36 registrations and their ordered
populations before fitting, runs one fit at a time, records phase peak RSS and
samples WSL memory, GPU and disk every ten seconds. Reuse the same matrix and
audit directories when resuming. Keep publication and scratch on the same
filesystem for the atomic rename.

```bash
(
  set -e
  step="$PWD/cpst/05_tsg_loss/05.1_loss_metrics"
  work=$(mktemp -d "$PWD/tmp/tsg-loss-formal.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
  capped_step04 python -u "$step/execution/matrix.py" \
    --output "$work/matrix" --logs "$work/matrix-audit"
  report_args=(--matrix "$work/matrix" --matrix-audit "$work/matrix-audit"
    --report "$work/report" --publication "$work/publication"
    --activity "$work/report-audit" --scratch "$work/scratch")
  capped_step04 python "$step/execution/report_formal.py" build \
    "${report_args[@]}" --code-revision "$(git rev-parse HEAD)"
  capped_step04 python "$step/execution/report_formal.py" verify "${report_args[@]}"
)
```

The matrix accepts `--training`, `--representations` and `--source`. The
report accepts `--training`, `--bandwidth-training`, `--previous-runs` and
`--previous-bandwidths` for relocated releases. Require exact agreement with
Step 04 for the three native `mse`, `alpha=1` anchors on epochs, objective
traces, predictions, Delta arrays, scores and eligibility counts. Verify the
complete ledger with the Step 04 metric order, frozen HVG-DEG 5k bandwidths
and equal-fold aggregation. [Formal receipt](log/330-formal.json) ·
[Fold score ledger](log/fold-scores.csv) · [Frozen bandwidths](log/bandwidths.csv)

## Replay a sealed fit

The formal registrations bind `997ca0f1ef5fd592b5faeb0979cb03cd5bf16f0c`
and paths under `stage2/`. Use the pinned runtime and original code in a
detached checkout. This example replays a scratch copy of native MSE, fold 0.
[Registration verification](../../util/training/verification.md) ·
[Historical publication and report verification](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/05_tsg_loss/05.1_loss_metrics/README.md)

```bash
(
  set -e
  root="$PWD"
  stage="$root/cpst"
  study="05_tsg_loss/05.1_loss_metrics"
  work=$(mktemp -d "$root/tmp/tsg-loss-pinned.XXXXXX")
  git worktree add --detach "$work/code" 997ca0f1ef5fd592b5faeb0979cb03cd5bf16f0c
  trap 'git worktree remove "$work/code"' EXIT
  mkdir "$work/run" "$work/scratch"
  sealed_run="$stage/$study/output/formal/matrix/runs/mse/alpha-1/fold-0"
  cp --reflink=auto --dereference -- "$sealed_run"/* "$work/run/"
  export TMPDIR="$work/scratch" PYTHONDONTWRITEBYTECODE=1
  capped_step04 python "$work/code/stage2/$study/execution/run_loss.py" replay \
    --output "$work/run" --training "$stage/03_data_train/03.7_tenv_denoise-5k/output" \
    --representations "$stage/util/output/representations" \
    --source "$root/pert_model/TranSiGen/source"
)
```

## Verify the current figures

For changes to the current figure data, rebuild Step 05.1 from sealed fits
through `build_comparison.py --step 05.1_loss_metrics --output <candidate>`
in the capped runtime. Use a fresh candidate under `tmp/` and pass
`--evidence-root` when the releases are materialized in another checkout.
[Accepted-fit rebuild](../README.md#rebuild-the-current-comparisons-from-accepted-fits)

Render the three published PNGs from their plotted CSVs and compare them byte
for byte. Use the [joint verifier](../README.md#verify-the-current-figures)
to check the receipts and two independent renders before replacing the figure
set. [Render receipt](log/figure-render.json)

```bash
(
  set -e
  step="$PWD/cpst/05_tsg_loss/05.1_loss_metrics"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-loss-figures.XXXXXX")
  TMPDIR="$work" MPLCONFIGDIR="$work/matplotlib" \
    python "$step/figure/script/render_comparison.py" --data "$step/figure/data" --output "$work/replayed"
  for name in learning-curve.png metrics.png metrics-paired.png; do
    cmp "$step/figure/$name" "$work/replayed/$name"
  done
)
```
