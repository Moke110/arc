# TranSiGen Delta loss comparison

[Spec 8](https://github.com/Moke110/arc/issues/324) crosses four Delta loss
forms (MSE, cosine, Pearson correlation and L1) with the weights 1/4, 1/3
and 1/2. It uses the [denoise-5k training release](../03_data_train/03.7_tenv_denoise-5k/README.md),
three frozen cold-drug folds, seed 42, and the same-plate control and TranSiGen
model from [Step 04](../04_tsg_axis/README.md). The held-out analysis measures
prediction conditional on a measured matched plate control.

A weight `alpha` gives the Delta term coefficient `3 * alpha` and each native
reconstruction MSE term coefficient `1.5 * (1 - alpha)`. The three native
KL terms each keep coefficient 0.1. The MSE × 1/3 coordinate calls the
native author objective. MSE and L1 sum residual penalties over conditions
and genes. Cosine and Pearson correlation sum per-condition dissimilarities
and multiply by the axis width. Pearson correlation first centers each Delta
vector over genes. Cosine and Pearson correlation clamp each L2 denominator
to at least `1e-8`. Both Delta vectors subtract the same measured control.
See [ADR 0088](../../docs/adr/0088-transigen-delta-loss-comparison.md).

Each fit uses its own objective to select the minimum validation loss with
30-epoch patience and a 1,000-epoch cap. The trace records loss per condition
and the six raw, unweighted objective components for both partitions.
Separate processes verify the chosen checkpoint, validation objective,
validation prediction and held-out prediction.

## Replay the bounded comparison

Activate the `arc` environment. Run from the materialized checkout root
using the [Step 04 capped runtime](../04_tsg_axis/README.md#prepare-the-capped-runtime).
Use the same environment and thread caps in every process.

```bash
set -e
STEP=stage2/05_tsg_loss
WORK=$(mktemp -d "$PWD/tmp/tsg-loss-smoke.XXXXXX")
runs=()
for loss in mse cos pcc l1; do
  run="$WORK/$loss"
  capped_step04 python "$STEP/execution/run_loss.py" register \
    --loss "$loss" --alpha 1/3 --fold 0 --smoke --output "$run"
  runs+=(--run "$run")
done
for loss in mse cos pcc l1; do
  for phase in train predict replay; do
    capped_step04 python "$STEP/execution/run_loss.py" "$phase" --output "$WORK/$loss"
  done
done
for phase in score verify; do
  capped_step04 python "$STEP/execution/evaluate_loss.py" "$phase" --smoke \
    "${runs[@]}" --output "$WORK/report"
done
```

Use `--training`, `--representations` and `--source` when an input is
materialized in another checkout. The evaluator accepts `--training` for
Denoise-DEG 5k and `--bandwidth-training` for HVG-DEG 5k. Input digests
and ordered conditions must still match the registered runs.

The [smoke receipt](log/325-smoke.json) records four two-epoch fits with
64 training, 24 validation and 24 test conditions each. It verifies the
native anchor, alternate objective calculations, exact model replay and
the nine Step 04 scores on the common held-out rows. Its original report
manifest is pinned to evaluator revision `3e4cbe32`.

## Run the formal matrix

The formal grid registers all 36 fits before the first training run. It
compares the full ordered train, validation and test populations of each arm
to the released folds. Runs use one CUDA fit at a time in a 10 GiB memory
and 8 GiB swap scope, with two CPU threads. The driver schedules WSL RAM,
GPU memory and disk samples on a monotonic 10-second cadence across phase
transitions. It records each phase's peak RSS.
An interrupted matrix can be resumed with the same command and candidate
paths. A completed phase keeps its registration and hashes.

```bash
set -e
STEP="$PWD/stage2/05_tsg_loss"
WORK=$(mktemp -d "$PWD/tmp/tsg-loss-formal.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
capped_step04 python -u "$STEP/execution/matrix.py" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
capped_step04 python "$STEP/execution/report_formal.py" build \
  --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
  --report "$WORK/report" --publication "$WORK/publication" \
  --activity "$WORK/report-audit" --scratch "$WORK/scratch" \
  --code-revision "$(git rev-parse HEAD)"
```

Use the same pinned inputs in each phase. If the run workspace has a separate
checkout, supply the materialized input paths as in the smoke instructions.
Use `--previous-runs` and `--previous-bandwidths` when the sealed Step 04
matched-control products live in another checkout.

## Nine held-out Delta metrics

[The four-panel figure](figure/delta-metrics.png) groups the Step 04 metrics
by loss form and uses three colors for the loss weights. Its source
[metric summary](figure/data/delta-metrics.csv) weights the three folds
equally. The [fold scores](log/fold-scores.csv) include condition
eligibility, and [run identities](log/runs.csv) bind the selected and
stopping epochs to the exact training and prediction evidence.
The [actual fold bandwidths](log/bandwidths.csv) match the Step 04 values.

The [formal receipt](log/326-formal.json) records the resource
[samples](log/326-resources.jsonl) and pins the local products under
`output/formal/matrix`, `output/formal/matrix-audit` and
`output/formal/report`. During the original fit, 35 intervals between
in-scope samples exceeded 12 seconds and the longest lasted 20.003 seconds.
All recorded samples met the resource limits. The MSE × 1/3 checkpoint predictions and objective
traces replay the three accepted Step 04 results exactly in the recorded
environment. Run `report_formal.py verify` with those local paths, a new
activity and scratch directory, and `--publication "$STEP"` to repeat
the end-to-end score and figure verification.

Render the metric figure from its source CSV alone.

```bash
mkdir -p "$PWD/tmp/tsg-loss-figure"
python stage2/05_tsg_loss/figure/script/render_metrics.py \
  --data stage2/05_tsg_loss/figure/data --output "$PWD/tmp/tsg-loss-figure"
```

## Learning curves

[The twelve-panel figure](figure/learning-curves.png) shows the training and
validation objective per condition for all 36 fits. Each fold keeps the same
color across panels. Dashed lines show training and solid lines show
validation. Every trace ends at its observed stopping epoch. The shared
logarithmic vertical scale keeps the early peaks and later changes visible.

The [epoch CSV](figure/data/epoch-objectives.csv) preserves all 4,817 accepted
epoch rows, including the objective components and selected and stopping
epochs. The [learning-curve receipt](log/327-curves.json) binds it to
the [formal receipt](log/326-formal.json) and the [run and checkpoint records](log/runs.csv).
The audit verifies contiguous epochs, checkpoint selection and stopping for
every run, then checks that an independent CSV-only render has the published
PNG digest.

```bash
set -e
STEP="$PWD/stage2/05_tsg_loss"
WORK=$(mktemp -d "$PWD/tmp/tsg-loss-curves.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
python "$STEP/figure/script/render_learning_curves.py" \
  --data "$STEP/figure/data" --output "$WORK"
python "$STEP/execution/audit_learning_curves.py" \
  --publication "$STEP" --scratch "$WORK/scratch" > "$WORK/327-curves.json"
cmp "$WORK/327-curves.json" "$STEP/log/327-curves.json"
```
