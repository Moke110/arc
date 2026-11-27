# TranSiGen additive Delta L1 sweep

[Spec 10](https://github.com/Moke110/carps/issues/331) sweeps nine additive
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
STEP="$PWD/stage2_CarPS-VAE/05_tsg_loss/05.2_mse_l1"
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

## Reproduce the original formal sweep

Use the same capped runtime. The matrix registers all 27 coordinates before
fitting, verifies the ordered conditions against the frozen folds and fits
one coordinate at a time. It records each phase's peak RSS and samples WSL
RAM, GPU and disk on a ten-second schedule. Resume with the same paths.
The matrix validates any sealed phase before proceeding.

```bash
set -e
STEP="$PWD/stage2_CarPS-VAE/05_tsg_loss/05.2_mse_l1"
WORK=$(mktemp -d "$PWD/tmp/tsg-mse-l1-formal.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
capped_step04 python -u "$STEP/execution/matrix_sweep.py" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
capped_step04 python "$STEP/execution/report_sweep.py" build \
  --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
  --report "$WORK/report" --publication "$WORK/publication" \
  --activity "$WORK/report-audit" --scratch "$WORK/scratch" \
  --code-revision "$(git rev-parse HEAD)"
```

Use `--training`, `--representations` and `--source` in both commands
when the same frozen inputs are materialized in another checkout. The
report command also accepts `--bandwidth-training` and `--previous-step`
to locate the frozen HVG-DEG 5k release and accepted Step 05.1 study.

## Verify the published study

The [formal receipt](log/333-formal.json) binds 27 independently replayed
fits, 243 [fold scores](log/fold-scores.csv), 81
[equal-fold scores](figure/data/delta-metrics.csv), the
[run ledger](log/runs.csv), and the [frozen bandwidths](log/bandwidths.csv).
All three zero-weight controls match Step 05.1 native MSE at alpha 1
exactly on the native objective trace, selected and stopping epochs,
validation predictions, test predictions, Delta arrays, fold scores and
eligibility counts. The [scoped resource log](log/333-resources.jsonl) contains 137
samples on a ten-second schedule. Twenty observed gaps exceeded 12 seconds,
with a maximum gap of 19.951 seconds. The
[independent coordinator log](log/333-coordinator-resources.jsonl) contains
137 GPU and WSL RAM samples across the same interval, with a maximum
monotonic gap of 10.147 seconds. The formal receipt binds both logs.

The [raw scores](figure/metrics.png) and
[signed score differences](figure/metrics-paired.png) show eight metrics, including wJAC-2000.
Bars show means, dots show individual [plotted fold scores](figure/data/metrics.csv), and error bars show sample SD.
The signed difference subtracts the matching zero-weight score within each
run. The [learning curves](figure/learning-curve.png) show one model per
alpha and a training and validation line per fold/seed run through its stopping
epoch. The [sealed objective CSV](figure/data/epoch-objectives.csv) retains the
objective components.

The [render receipt](log/figure-render.json) binds the figures to the
accepted release.
[Current figure verification](../README.md#verify-the-current-figures)
compares the refreshed PNGs with two independent CSV-only renders.

Render all three PNGs from the published CSVs alone.

```bash
STEP="$PWD/stage2_CarPS-VAE/05_tsg_loss/05.2_mse_l1"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-loss-figures.XXXXXX")
export TMPDIR="$WORK" MPLCONFIGDIR="$WORK"
python "$STEP/figure/script/render_comparison.py" --data "$STEP/figure/data" --output "$WORK"
```
