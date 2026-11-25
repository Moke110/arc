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
mkdir -p "$PWD/tmp"
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

## Run the formal sweep

Use the capped runtime for a single CUDA fit at a time. The matrix registers
all 27 coordinates before fitting and checks sealed phase files before
resuming on the same paths. Its audit samples WSL RAM and GPU usage on a
ten-second schedule and records peak RSS for each process.

Allow at least 4 GiB free on the volume used for formal work. The matrix
takes about 3 GiB and the runner requires 1 GiB of free space during fitting.
To use a different volume, set `TSG_RECONSTRUCTION_WORK_ROOT` to an existing
directory or symbolic link under this checkout's `tmp/` backed by that volume.
Keep the scratch, report and publication directories on the same filesystem.
The report builder promotes staged directories with filesystem renames.

```bash
set -e
STEP="$PWD/stage2/05_tsg_loss/05.3_mse-l1_mse-x"
mkdir -p "$PWD/tmp"
WORK_ROOT="${TSG_RECONSTRUCTION_WORK_ROOT:-$PWD/tmp}"
WORK=$(mktemp -d "$WORK_ROOT/tsg-reconstruction-formal.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
REVISION=$(git rev-parse HEAD)
capped_step04 python -u "$STEP/execution/matrix_reconstruction.py" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
capped_step04 python "$STEP/execution/report_reconstruction.py" build \
  --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
  --report "$WORK/report" --publication "$WORK/publication" \
  --activity "$WORK/report-audit" --scratch "$WORK/scratch" \
  --code-revision "$REVISION"
```

The report can use `--previous-step` to locate the sealed Step 05.2 release
in another checkout. It also accepts the runner input options and
`--bandwidth-training`. Use the same report command with `verify` and omit
`--code-revision` to repeat the report acceptance checks.

## Verify the published sweep

The accepted release contains 27 independently replayed runs, 243 held-out
fold scores and 81 equal-fold aggregates. The acceptance receipt binds
the matrix, run ledger, ordered fold membership, frozen bandwidths and
resource records. [Acceptance receipt](log/336-formal.json) ·
[Run ledger](log/runs.csv) · [Fold scores](log/fold-scores.csv) ·
[Bandwidths](log/bandwidths.csv)

The coordinator recorded 125 GPU and WSL RAM samples with a maximum
monotonic gap of 10.002 seconds. The scoped runner recorded 124 samples,
with a maximum gap of 13.232 seconds, and recorded peak RSS for each phase.
[Coordinator samples](log/336-coordinator-resources.jsonl) ·
[Scoped runner samples](log/336-resources.jsonl)

All three alpha-one fits reproduced Step 05.2 at Delta L1 weight 1.75
exactly on the train and validation objective and component traces,
selected and stopping epochs, validation and test predictions, Delta
arrays, fold scores and eligibility counts. The nine alpha-one signed
differences are zero.

The score and signed-difference plots show the same nine held-out metrics.
Bars show means, dots show the individual [fold scores](figure/data/fold-scores.csv),
and error bars show sample SD. Differences subtract the matching alpha-one
score within each run. Learning-curve panels show one model per reconstruction
weight and one training and validation line per fold/seed run.
[Score figure](figure/delta-metrics.png) ·
[Signed differences](figure/delta-metrics-improvement.png) ·
[Learning curves](figure/learning-curves.png) ·
[Metric CSV](figure/data/delta-metrics.csv) ·
[Learning-curve CSV](figure/data/epoch-objectives.csv)

The [current figure verification](../README.md#verify-the-current-figures)
compares the refreshed PNGs with two independent CSV-only renders.
