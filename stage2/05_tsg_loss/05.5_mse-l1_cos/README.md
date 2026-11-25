# TranSiGen Delta cosine sweep

[Spec 12](https://github.com/Moke110/arc/issues/337) adds a bounded cosine similarity
loss on predicted and measured Delta from the same matched-plate DMSO
control. The sweep fixes the Step 05.3 shared reconstruction weight at 1
and Delta L1 weight at 1.75. The added loss is the batch mean of half of
one minus similarity, with each similarity clipped to [−1, 1].
The six weights are 0, 0.1, 0.2, 0.3, 0.4 and 0.5.
[ADR 0092](../../../docs/adr/0092-transigen-bounded-delta-similarity-sweeps.md)

## Reproduce the smoke

Use the `arc` environment and [capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime)
from the materialized checkout. Monitor WSL RAM and GPU usage every ten
seconds. The [smoke receipt](log/338-smoke.json) covers two fold-zero
fits with two epochs each, exact prediction replay, and two identical
renders from the report CSVs alone.

```bash
set -e
metric=cos
LOSS="$PWD/stage2/05_tsg_loss"
SHARED="$LOSS/similarity"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-$metric-smoke.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
runs=()
for alpha in 0 0.5; do
  run="$WORK/alpha-$alpha"
  capped_step04 python "$SHARED/run_similarity.py" register \
    --metric "$metric" --alpha "$alpha" --fold 0 --smoke --output "$run"
  runs+=(--run "$run")
done
for alpha in 0 0.5; do
  for phase in train predict replay; do
    capped_step04 python "$SHARED/run_similarity.py" "$phase" \
      --metric "$metric" --output "$WORK/alpha-$alpha"
  done
done
for phase in score verify; do
  capped_step04 python "$SHARED/evaluate_similarity.py" "$phase" \
    --metric "$metric" --smoke "${runs[@]}" --output "$WORK/report" --scratch "$WORK/scratch"
done
mkdir "$WORK/csv-only"
cp "$WORK/report/"{delta-metrics,epoch-objectives}.csv "$WORK/csv-only/"
for attempt in first second; do
  capped_step04 python "$LOSS/05.5_mse-l1_cos/figure/script/render.py" --smoke \
    --data "$WORK/csv-only" --output "$WORK/$attempt"
done
for figure in delta-metrics.png delta-metrics-improvement.png learning-curves.png; do
  cmp "$WORK/first/$figure" "$WORK/second/$figure"
done
```

Use `--training`, `--representations` and `--source` to locate the
pinned runner inputs in another checkout. The evaluator accepts
`--training` and `--bandwidth-training`.

## Verify the published sweep

The release contains 18 fresh CUDA fits with exact independent replay, 162
held-out fold scores and 54 equal-fold aggregates. Its three zero-weight
fits agree exactly with the corresponding Step 05.3 unit-weight references
on the objective and component traces, selected and stopping epochs,
validation and test predictions, Delta arrays, nine scores and eligibility
counts. The acceptance receipt binds the shared matrix and registered
inputs, the exact reference comparisons, and the scoring and resource
evidence. [Acceptance receipt](log/339-formal.json) · [Run ledger](log/runs.csv) ·
[Fold scores](log/fold-scores.csv) · [Frozen bandwidths](log/bandwidths.csv)

Score bars show means, dots show the individual
[fold scores](figure/data/fold-scores.csv), and error bars show sample SD.
Differences subtract the matching fresh zero-weight score within each run.
The learning curves show one model per weight and one training and validation
line per fold/seed run. Dots mark the selected checkpoints.
[Scores](figure/delta-metrics.png) ·
[Signed differences](figure/delta-metrics-improvement.png) ·
[Learning curves](figure/learning-curves.png) ·
[Score CSV](figure/data/delta-metrics.csv) ·
[Curve CSV](figure/data/epoch-objectives.csv)

The independent coordinator recorded 172 samples with a maximum observed
gap of 10.002 seconds. The capped fit scope recorded 170 samples with a
maximum gap of 14.421 seconds and peak RAM of 1.84 GiB.
[Coordinator samples](log/339-coordinator-resources.jsonl) ·
[Scope samples](log/339-resources.jsonl)

The [current figure verification](../README.md#verify-the-current-figures)
compares the published PNGs with two independent CSV-only renders.
