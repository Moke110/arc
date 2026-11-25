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

## Run the formal sweep

Use the same capped runtime. The matrix registers all 27 coordinates before
fitting, verifies the ordered conditions against the frozen folds and fits
one coordinate at a time. It records each phase's peak RSS and samples WSL
RAM, GPU and disk on a ten-second schedule. Resume with the same paths.
The matrix validates any sealed phase before proceeding.

```bash
set -e
STEP="$PWD/stage2/05_tsg_loss/05.2_mse_l1"
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
eligibility counts. The [resource samples](log/333-resources.jsonl) record
the monitored execution.

The [raw scores](figure/delta-metrics.png) and
[signed score differences](figure/delta-metrics-improvement.png) show one
panel per metric. The difference is the score at the plotted alpha minus
the zero-weight score. Blue backgrounds identify metrics where higher
values are better. Red backgrounds identify metrics where lower values
are better. The [learning curves](figure/learning-curves.png) show the
training and validation objectives through the stopping epoch, with one
panel per alpha and one color per fold. The
[epoch CSV](figure/data/epoch-objectives.csv) retains the raw objective
components.

To verify the active result, run the formal report command above with
`verify`, the `matrix`, `matrix-audit` and `report` directories under
`$STEP/output/formal`, and `--publication "$STEP"`. Use fresh scratch and
activity directories and omit `--code-revision`.

Render all three PNGs from the published CSVs alone.

```bash
STEP="$PWD/stage2/05_tsg_loss/05.2_mse_l1"
WORK=$(mktemp -d "$PWD/tmp/tsg-mse-l1-figures.XXXXXX")
export TMPDIR="$WORK" MPLCONFIGDIR="$WORK"
for renderer in render_metrics render_learning_curves; do
  python "$STEP/figure/script/$renderer.py" --data "$STEP/figure/data" --output "$WORK"
done
```
