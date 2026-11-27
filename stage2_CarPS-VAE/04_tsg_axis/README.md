# TranSiGen prediction-axis comparison

The numbered commands use the [shared TranSiGen training and evaluation APIs](../util/README.md#transigen-training-and-evaluation).

[Spec 6](https://github.com/Moke110/carps/issues/317) compares three HVG-DEG
and three Denoise-DEG axes on the same tenv plate–parent–dose conditions.
Each fit starts with random weights at seed 42 and uses the original
TranSiGen internal widths, frozen 2,304-dimensional KPGT and the registered
latent numeric dose transform. The selected prediction axis sets the gene
interface width. Fold membership comes from the frozen cold-drug folds in
the [Step 03 tenv training release](../03_data_train/README.md).

Training, validation and test each use the released same-plate DMSO expression
matched to the condition. Training and validation use the author-published
TranSiGen MSE + KLD objective with beta 0.1. Validation sums the objective
over batches in evaluation mode with a fixed seed and divides by the number
of validation conditions. The minimum selects the checkpoint. Training stops
after 30 consecutive epochs without improvement or at 1,000 epochs.

At the fixed prediction boundary, `prediction.npy` stores absolute treated
expression. Predicted and measured Delta both subtract the same matched plate
control. All three matrices follow `conditions.csv` and `genes.csv`. Test
treated expression is used for scoring after prediction. The control and
validation protocol is recorded in [ADR 0086](../../docs/adr/0086-tenv-axis-matched-plate-control.md).

## Gene axis selection

Select **Denoise-DEG 5k (`denoise-5k`)** as the Stage 2 prediction axis.
The decision prioritizes held-out perturbation Delta prediction with a measured
same-plate control. Across the six axes and three frozen cold-drug folds,
Denoise-DEG 5k gave the highest mean Δ PCC (0.4066) and the lowest mean
Δ RMSE (0.2596) and Δ MMD (0.1793). HVG-DEG variants scored higher on signed
DEG F1 for each of the four disease references. See the
[formal results](#formal-results) for the full comparison.

## Prepare the capped runtime

Activate the `arc` environment and run from the repository root. Use one CUDA
job at a time. The launcher checks for at least 12 GiB of available WSL RAM
before each command and sets the measured run's limits of 10 GiB of memory,
8 GiB of swap and two CPU threads. Monitor available RAM and GPU usage every
10 seconds while fitting.

```bash
capped_step04() {
  if ! LC_ALL=C free -g | awk '
    { print }
    $1 == "Mem:" { available = $7 }
    END { exit available < 12 }
  '; then
    printf '%s\n' 'Step 04 requires at least 12 GiB of available RAM.' >&2
    return 1
  fi
  systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- \
    env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 "$@"
}
```

## Run the formal matrix

The matrix driver checks the complete ordered condition and fold rosters
across all six releases. It registers all 18 coordinates before training,
then runs each fit, prediction and exact replay in separate processes.
For the formal dataset, use the full cold-drug folds and seed 42.

```bash
STEP=stage2_CarPS-VAE/04_tsg_axis
WORK=$(mktemp -d "$PWD/tmp/tsg-axis-matched-control.XXXXXX")
capped_step04 python "$STEP/execution/matrix.py" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit/logs" &&
capped_step04 python "$STEP/execution/report_formal.py" build \
  --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
  --report "$WORK/report" --figures "$WORK/figures" \
  --logs "$WORK/logs"
```

The accepted [formal receipt](log/319-formal.json) pins all 18 run
manifests, the report and the original two figures. The materialized
paths are `output/formal-matched-control` for checkpoints, registrations
and predictions, `log/319-runs` for resource and process evidence, and
`output/report-matched-control` for condition and fold metrics, eligibility,
bandwidths and the shared disease reference. This sealed report keeps its
original axis IDs and digests. The current plotted CSVs in `figure/data/` use canonical IDs.

Export the six axis-bearing condition, fold and summary CSVs with canonical
IDs. The export verifies the sealed report's manifest and output digests
against the formal receipt and copies each row's other fields byte for byte.

```bash
STEP=stage2_CarPS-VAE/04_tsg_axis
WORK=$(mktemp -d "$PWD/tmp/tsg-axis-export.XXXXXX")
python "$STEP/execution/export_canonical.py" --output "$WORK/report"
```

Verify the sealed scientific report with its pinned evaluator and axis IDs.
The temporary Step 03 aliases point to the same releases.

```bash
STEP=stage2_CarPS-VAE/04_tsg_axis
VERIFY=$(mktemp -d "$PWD/tmp/tsg-axis-verify.XXXXXX")
git worktree add --detach "$VERIFY/runtime" 3f4cbb65
mkdir "$VERIFY/legacy-step03"
for release in "$PWD"/stage2_CarPS-VAE/03_data_train/03.*_tenv_*; do
  name=${release##*/}
  ln -s "$release" "$VERIFY/legacy-step03/${name%-*}-deg-${name##*-}"
done
runs=()
for run in "$STEP"/output/formal-matched-control/runs/*/fold-*; do
  runs+=(--run "$run")
done
capped_step04 python "$VERIFY/runtime/stage2/04_tsg_axis/execution/evaluate.py" verify "${runs[@]}" \
  --output "$STEP/output/report-matched-control" \
  --step03 "$VERIFY/legacy-step03" --raw "$PWD/stage2_CarPS-VAE/01_data_raw/output"
git worktree remove "$VERIFY/runtime"
```

## Register and replay one fit

```bash
STEP=stage2_CarPS-VAE/04_tsg_axis
RUN=$PWD/tmp/tsg-axis/hvg-5k-fold-0
capped_step04 python "$STEP/execution/run.py" register --axis hvg-5k --fold 0 \
  --training stage2_CarPS-VAE/03_data_train/03.1_tenv_hvg-5k/output --device cuda --output "$RUN" &&
capped_step04 python "$STEP/execution/run.py" train --output "$RUN" &&
capped_step04 python "$STEP/execution/run.py" predict --output "$RUN" &&
capped_step04 python "$STEP/execution/run.py" replay --output "$RUN"
```

Pass `--representations` and `--source` when frozen KPGT or immutable native
TranSiGen is materialized in a different checkout. Keep the same CPU thread
caps for each phase. Registration binds dataset and gene order, parent
and condition partitions, feature and native code digests and the runtime.
Prediction exactly replays the validation loss and validation predictions
before producing held-out predictions.
A fresh replay verifies all scores' source arrays and row identities against
the recorded values and publishes the run manifest. Repeating replay checks
the sealed run again without changing it. Pass the corresponding path flags
when an input has moved. Its registered digests still apply.

Pass `--smoke` only for a bounded two-epoch fit with 64 training, 24
validation and 24 test conditions. Its manifest marks it as smoke and it is
excluded from formal results. The [smoke receipt](log/318-matched-control-smoke.json)
records the matched-control 5k fit, scoring and figure replay.

## Score and replay an axis comparison

Pass one `--run` per sealed fit. The formal report requires all 18 axis–fold
coordinates. Add `--smoke` to score a single replayed HVG-DEG 5k smoke fit.
The evaluator checks all registered partitions and prediction arrays against
the Stage 03 releases. Use `--step03` and `--raw` for materialized inputs in
another checkout.

```bash
STEP=stage2_CarPS-VAE/04_tsg_axis
REPORT=$PWD/tmp/tsg-axis/report
python "$STEP/execution/evaluate.py" score --smoke --run "$RUN" --output "$REPORT"
python "$STEP/execution/evaluate.py" verify --smoke --run "$RUN" --output "$REPORT"
```

The report orders the metrics as JAC_top200, JAC_pos100, JAC_neg100, Δ PCC,
Δ Spearman, Δ RMSE, Δ E-distance, Δ MMD and Δ R². Jaccard counts strictly
positive and strictly negative genes separately and breaks ties by frozen gene
index. E-distance and RBF MMD use Euclidean distance divided by the square
root of axis width and include the empirical diagonal terms. Each fold uses
one kernel bandwidth for all six axes. It is the median over distinct
condition pairs in the full HVG-DEG 5k training measured Delta. The distance
and kernel bandwidth evidence is in `bandwidths.csv` and its condition
membership is in `fold-memberships.csv`.

For each eligible disease, signed DEG F1 compares the full model-axis
predictions against nonzero directions in its DEG1000 list, measured from
raw pooled counts and matched plate controls on the full common CM axis.
Opposite signs count as one false positive and one false negative. An off-axis
reference gene counts as a false negative. Two empty sets produce NA. The
same `reference/` files supply all models' disease targets.

The report retains per-condition, per-fold and summary CSVs with eligibility counts.
The report weights defined conditions equally within a fold and defined folds
equally within each axis. The current comparison also includes wJAC-2000, a
weighted average of the signed Jaccard cutoffs 200, 500, 1,000 and 2,000 per
condition. Its fold scores use the same condition aggregation. These are
retrospective results on published fixed axes built with tenv evidence and
measure prediction conditional on a measured matched control.

## Formal results

All 18 fits passed exact replay of validation loss, validation predictions
and held-out predictions. Training stopped after 127 to 243 epochs, with
30 epochs without validation loss improvement in every fit. All six axes
evaluated the same 462, 459 and 462 conditions in folds 0, 1 and 2.
The fold RBF bandwidths were 0.430414, 0.425762 and 0.439392. Every metric
and every disease score was defined in all three folds.

| Delta metric | `hvg-3k` | `hvg-4k` | `hvg-5k` | `denoise-3k` | `denoise-4k` | `denoise-5k` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| JAC_top200 | 0.1453 | 0.1428 | 0.1416 | 0.1454 | 0.1427 | 0.1409 |
| JAC_pos100 | 0.1375 | 0.1319 | 0.1258 | 0.1283 | 0.1227 | 0.1134 |
| JAC_neg100 | 0.1644 | 0.1650 | 0.1695 | 0.1683 | 0.1700 | 0.1766 |
| Δ PCC | 0.3928 | 0.3967 | 0.4043 | 0.4023 | 0.4003 | 0.4066 |
| Δ Spearman | 0.3896 | 0.3989 | 0.4045 | 0.3906 | 0.3888 | 0.3961 |
| Δ RMSE | 0.3586 | 0.3240 | 0.2957 | 0.3107 | 0.2893 | 0.2596 |
| Δ E-distance | 0.0370 | 0.0394 | 0.0350 | 0.0320 | 0.0315 | 0.0319 |
| Δ MMD | 0.2253 | 0.2208 | 0.1994 | 0.1955 | 0.1873 | 0.1793 |
| Δ R² | -0.1953 | -0.2407 | -0.1842 | -0.1993 | -0.2275 | -0.1587 |

| Signed disease DEG F1 | `hvg-3k` | `hvg-4k` | `hvg-5k` | `denoise-3k` | `denoise-4k` | `denoise-5k` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nHF | 0.1603 | 0.1639 | 0.1905 | 0.1197 | 0.1063 | 0.0969 |
| ICM | 0.1083 | 0.1258 | 0.1734 | 0.1211 | 0.1057 | 0.0961 |
| genetic cardiomyopathy | 0.2136 | 0.2469 | 0.2193 | 0.1208 | 0.1057 | 0.0950 |
| TOX | 0.3026 | 0.2587 | 0.2162 | 0.1302 | 0.1148 | 0.1036 |

See the [learning curves](figure/learning-curve.png), the
[held-out scores](figure/metrics.png) and the [paired differences](figure/metrics-paired.png).
The paired figure subtracts Denoise-DEG 5k from each other axis in the same
fold and seed. Every metric panel includes zero, the defined fold scores,
their mean and sample SD. The
[previous publication](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/04_tsg_axis/README.md#formal-results)
contains the disease F1 figure and historical rendering receipt.

Rebuild the comparison directly from the sealed Step 04 fits and report.
The build checks all 18 registrations, rehashes the matched-control Delta
arrays and epoch traces, and checks JAC-top200 against the sealed
per-condition and fold scores.

```bash
STEP=stage2_CarPS-VAE/04_tsg_axis
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-axis-figures.XXXXXX")
export TMPDIR="$WORK" MPLCONFIGDIR="$WORK/mpl"
python "$STEP/execution/build_comparison.py" --output "$WORK/publication"
```

Reproduce the three PNGs using only the plotted CSVs.

```bash
python "$STEP/figure/script/render.py" --data "$STEP/figure/data" \
  --output "$WORK/replayed"
```

The [current rendering receipt](log/figure-render.json) pins the
[trace CSV](figure/data/learning-curve.csv), [score CSV](figure/data/metrics.csv),
[paired CSV](figure/data/metrics-paired.csv) and all three reproducible PNGs.
The [canonical export receipt](log/322-canonical-axis-report.json) records the
historical six ID mappings and sealed table hashes.

## Historical comparisons

The [previous routing](https://github.com/Moke110/carps/blob/59b188ec0fe6386f2c916e062eaad310c20f00f8/stage2/04_tsg_axis/README.md)
records the original inputs and pinned replay commands for the single-dose and
axis-wide-control comparisons and the subsequent dose and width studies.
Use the [original tracked evidence](../README.md#historical-evidence-and-active-input-checks)
for the earlier comparisons.