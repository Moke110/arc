# TranSiGen Delta loss metrics

[Spec 9](https://github.com/Moke110/arc/issues/328) compares four Delta loss
metrics (MSE, cosine, Pearson correlation and L1) at coefficients 0.75, 1 and
1.25. The study uses the [denoise-5k training release](../../03_data_train/03.7_tenv_denoise-5k/README.md),
three frozen cold-drug folds, seed 42 and the same-plate control and TranSiGen
model from [Step 04](../../04_tsg_axis/README.md). The held-out analysis measures
prediction conditional on a measured matched plate control.

The objective is `0.1 * K + MSE_X1 + MSE_X2 + alpha * metric_Delta`, where
`K` is the sum of the three native KL terms. Both reconstruction MSE
coefficients are 1. MSE at alpha 1 calls the native author objective.
MSE and L1 sum over conditions and genes. Cosine and Pearson correlation sum
per-condition dissimilarities and multiply by the gene axis width. Pearson
correlation first centers each Delta vector over genes. Both angular losses
clamp the L2 denominators to at least `1e-8`. Both Delta vectors subtract the
same measured control. See [ADR 0089](../../../docs/adr/0089-transigen-delta-metric-weight.md).

Each fit selects the minimum validation loss on its registered objective
with 30-epoch patience and a 1,000-epoch cap. The trace records training and
validation losses per condition and their six raw objective components.
Separate processes verify the chosen checkpoint, validation objective,
validation prediction and held-out prediction.

## Run a bounded comparison

Activate the `arc` environment and define the
[Step 04 capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime)
from the materialized checkout root.

```bash
set -e
STEP="$PWD/stage2/05_tsg_loss/05.1_loss_metrics"
WORK=$(mktemp -d "$PWD/tmp/tsg-loss-smoke.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
runs=()
for loss in mse cos pcc l1; do
  run="$WORK/$loss"
  capped_step04 python "$STEP/execution/run_loss.py" register \
    --loss "$loss" --alpha 1 --fold 0 --smoke --output "$run"
  runs+=(--run "$run")
done
for loss in mse cos pcc l1; do
  for phase in train predict replay; do
    capped_step04 python "$STEP/execution/run_loss.py" "$phase" --output "$WORK/$loss"
  done
done
for phase in score verify; do
  capped_step04 python "$STEP/execution/evaluate_loss.py" "$phase" --smoke \
    "${runs[@]}" --output "$WORK/report" --scratch "$WORK/scratch"
done
```

The [smoke receipt](log/329-smoke.json) describes the four two-epoch fits.
Use `--training`, `--representations` and `--source` for materialized inputs
in another checkout. The evaluator also accepts `--bandwidth-training` for
HVG-DEG 5k. Input digests and ordered conditions must match the registered runs.

## Run the formal matrix

The formal driver registers all 36 fits before the first training run and
checks their ordered train, validation and test populations against the frozen
folds. It runs one CUDA fit at a time in the capped scope, records each phase's
peak RSS and samples WSL memory, GPU and disk every ten seconds. Resume an
interrupted run with the same paths to keep completed phases and their hashes.

```bash
set -e
STEP="$PWD/stage2/05_tsg_loss/05.1_loss_metrics"
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

Use `--previous-runs` and `--previous-bandwidths` when the accepted Step 04
matched-control products are materialized in another checkout. Keep candidate
publication and scratch on the same filesystem so publication can use an
atomic rename.

## Held-out scores and learning curves

[Raw scores](figure/delta-metrics.png) and
[signed score differences](figure/delta-metrics-improvement.png) each have
one panel per Step 04 held-out metric. The [metric CSV](figure/data/delta-metrics.csv)
includes the equal-fold score, the corresponding MSE at alpha 1 score and
their signed difference. Blue backgrounds indicate metrics where higher
scores are better. Red backgrounds indicate metrics where lower scores are
better. [Fold scores](log/fold-scores.csv) preserve eligibility counts and
[bandwidths](log/bandwidths.csv) record the frozen Step 04 values.

[Learning curves](figure/learning-curves.png) show the training objective
(dashed) and validation objective (solid) through the stopping epoch, with
one panel per loss and coefficient and one color per fold.
The [epoch CSV](figure/data/epoch-objectives.csv) holds the observed training
and validation values, the [run ledger](log/runs.csv) records the selected
checkpoint for each fit, and the [curve receipt](log/330-curves.json) binds
them to the published figure.

The [formal receipt](log/330-formal.json) binds all three native Step 04
anchors, the [resource samples](log/330-resources.jsonl), three PNGs and the
local products under `output/formal/{matrix,matrix-audit,report}`.
To verify the active result from the checkout root, run the formal report
command above with `verify`, those local matrix, audit and report paths,
and `--publication "$STEP"`. Use fresh scratch and activity paths and omit
`--code-revision`.

Render all three PNGs using only the published CSVs.

```bash
STEP="$PWD/stage2/05_tsg_loss/05.1_loss_metrics"
WORK=$(mktemp -d "$PWD/tmp/tsg-loss-figures.XXXXXX")
export TMPDIR="$WORK" MPLCONFIGDIR="$WORK"
python "$STEP/figure/script/render_metrics.py" --data "$STEP/figure/data" --output "$WORK"
python "$STEP/figure/script/render_learning_curves.py" --data "$STEP/figure/data" --output "$WORK"
```
