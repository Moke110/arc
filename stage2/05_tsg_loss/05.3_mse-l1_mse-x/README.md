# TranSiGen reconstruction MSE sweep

[Spec 11](https://github.com/Moke110/arc/issues/334) varies the shared
expression reconstruction MSE weight on the [denoise-5k training release](../../03_data_train/03.7_tenv_denoise-5k/README.md)
while fixing the Delta L1 weight at 1.75. The objective is
`0.1 * K + MSE_Delta + 1.75 * L1_Delta + alpha * (MSE_X1 + MSE_X2)`.
Both Delta vectors use the measured same-plate DMSO control. Each run
selects the minimum registered validation objective with 30-epoch
patience and a 1,000-epoch cap. At alpha one, the objective uses the exact
[Step 05.2](../05.2_mse_l1/README.md) arithmetic. Score differences subtract
the fresh alpha-one score from the score at each reconstruction weight.
See [ADR 0091](../../../docs/adr/0091-transigen-reconstruction-mse-sweep.md).

## Run the bounded comparison

Activate the `arc` environment and define the
[Step 04 capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime)
from the materialized checkout root. Monitor WSL RAM and GPU usage every
ten seconds. The [smoke receipt](log/335-smoke.json) records the bounded
two-epoch comparison on fold 0.

```bash
set -e
STEP="$PWD/stage2/05_tsg_loss/05.3_mse-l1_mse-x"
WORK=$(mktemp -d "$PWD/tmp/tsg-reconstruction-smoke.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
runs=()
for alpha in 0 1; do
  run="$WORK/alpha-$alpha"
  capped_step04 python "$STEP/execution/run_reconstruction.py" register \
    --alpha "$alpha" --fold 0 --smoke --output "$run"
  runs+=(--run "$run")
done
for alpha in 0 1; do
  for phase in train predict replay; do
    capped_step04 python "$STEP/execution/run_reconstruction.py" "$phase" \
      --output "$WORK/alpha-$alpha"
  done
done
for phase in score verify; do
  capped_step04 python "$STEP/execution/evaluate_reconstruction.py" "$phase" \
    --smoke "${runs[@]}" --output "$WORK/report" --scratch "$WORK/scratch"
done
mkdir "$WORK/csv-only"
cp "$WORK/report/"{delta-metrics,epoch-objectives}.csv "$WORK/csv-only/"
for attempt in first second; do
  for renderer in render_metrics render_learning_curves; do
    capped_step04 python "$STEP/figure/script/$renderer.py" --smoke \
      --data "$WORK/csv-only" --output "$WORK/$attempt"
  done
done
for figure in delta-metrics.png delta-metrics-improvement.png learning-curves.png; do
  cmp "$WORK/first/$figure" "$WORK/second/$figure"
done
```

Use `--training`, `--representations` and `--source` to locate the
pinned runner inputs in another checkout. The evaluator accepts
`--training` and `--bandwidth-training`.
