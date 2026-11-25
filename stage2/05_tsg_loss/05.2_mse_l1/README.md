# TranSiGen additive Delta L1 sweep

[Spec 10](https://github.com/Moke110/arc/issues/331) sweeps nine additive
Delta L1 weights from 0 to 2 in 0.25 increments across the three frozen
cold-drug folds on the [denoise-5k training release](../../03_data_train/03.7_tenv_denoise-5k/README.md).
Both Delta vectors use the measured same-plate DMSO control. The objective is
`0.1 * K + MSE_X1 + MSE_X2 + MSE_Delta + alpha * L1_Delta`, with the
native arithmetic preserved at weight zero and both raw Delta terms
recorded as per-condition batch sums. See [ADR 0090](../../../docs/adr/0090-transigen-mse-l1-sweep.md).

Each coordinate selects the minimum validation objective using 30-epoch
patience and a 1,000-epoch cap. Held-out scores use the Step 05.1 metric
definitions and frozen HVG-DEG 5k bandwidths. Score differences use the
fresh zero-weight runs as the reference.

## Run a bounded comparison

Activate the `arc` environment and define the
[Step 04 capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime)
from the materialized checkout root. Monitor WSL RAM and GPU usage every
ten seconds. The [smoke receipt](log/332-smoke.json) records two CUDA fits
on fold 0 and their independent replay.

```bash
set -e
STEP="$PWD/stage2/05_tsg_loss/05.2_mse_l1"
WORK=$(mktemp -d "$PWD/tmp/tsg-mse-l1-smoke.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
runs=()
for alpha in 0 1; do
  run="$WORK/alpha-$alpha"
  capped_step04 python "$STEP/execution/run_sweep.py" register \
    --alpha "$alpha" --fold 0 --smoke --output "$run"
  runs+=(--run "$run")
done
for alpha in 0 1; do
  for phase in train predict replay; do
    capped_step04 python "$STEP/execution/run_sweep.py" "$phase" \
      --output "$WORK/alpha-$alpha"
  done
done
for phase in score verify; do
  capped_step04 python "$STEP/execution/evaluate_sweep.py" "$phase" --smoke \
    "${runs[@]}" --output "$WORK/report" --scratch "$WORK/scratch"
done
for attempt in first second; do
  for renderer in render_metrics render_learning_curves; do
    capped_step04 python "$STEP/figure/script/$renderer.py" --smoke \
      --data "$WORK/report" --output "$WORK/$attempt"
  done
done
for figure in delta-metrics.png delta-metrics-improvement.png learning-curves.png; do
  cmp "$WORK/first/$figure" "$WORK/second/$figure"
done
```

Use `--training`, `--representations` and `--source` for the same
materialized runner inputs in another checkout. The evaluator accepts
`--training` and `--bandwidth-training`.
