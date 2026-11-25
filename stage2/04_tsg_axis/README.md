# TranSiGen prediction-axis comparison

[Spec 6](https://github.com/Moke110/arc/issues/317) compares three HVG-DEG
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
STEP=stage2/04_tsg_axis
WORK=$(mktemp -d "$PWD/tmp/tsg-axis-matched-control.XXXXXX")
capped_step04 python "$STEP/execution/matrix.py" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit/logs" &&
capped_step04 python "$STEP/execution/report_formal.py" build \
  --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
  --report "$WORK/report" --figures "$WORK/figures" \
  --logs "$WORK/logs"
```

The accepted [formal receipt](log/319-formal.json) pins all 18 run
manifests, the report and its two CSV-reproducible PNGs. The materialized
paths are `output/formal-matched-control` for checkpoints, registrations
and predictions, `log/319-runs` for resource and process evidence, and
`output/report-matched-control` for condition and fold metrics, eligibility,
bandwidths and the shared disease reference. The two summary CSVs are in
`figure/data/`.

```bash
STEP=stage2/04_tsg_axis
VERIFY=$(mktemp -d "$PWD/tmp/tsg-axis-verify.XXXXXX")
capped_step04 python "$STEP/execution/report_formal.py" verify \
  --matrix "$STEP/output/formal-matched-control" --matrix-audit "$STEP/log/319-runs" \
  --report "$STEP/output/report-matched-control" --figures "$STEP/figure" \
  --logs "$VERIFY"
```

## Register and replay one fit

```bash
STEP=stage2/04_tsg_axis
RUN=$PWD/tmp/tsg-axis/hvg-deg-5k-fold-0
capped_step04 python "$STEP/execution/run.py" register --axis hvg-deg-5k --fold 0 \
  --training stage2/03_data_train/03.1_tenv_hvg-deg-5k/output --device cuda --output "$RUN" &&
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
STEP=stage2/04_tsg_axis
REPORT=$PWD/tmp/tsg-axis/report
python "$STEP/execution/evaluate.py" score --smoke --run "$RUN" --output "$REPORT"
python "$STEP/execution/evaluate.py" verify --smoke --run "$RUN" --output "$REPORT"
python "$STEP/figure/script/render.py" --data "$REPORT/figure/data" \
  --output "$PWD/tmp/tsg-axis/figures"
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

The report retains per-condition and per-fold CSVs, eligibility counts and
the two aggregate CSVs consumed by the figures. It weights defined conditions
equally within a fold and defined folds equally within each axis. These are
retrospective results on published fixed axes built with tenv evidence and
measure prediction conditional on a measured matched control.

## Formal results

All 18 fits passed exact replay of validation loss, validation predictions
and held-out predictions. Training stopped after 127 to 243 epochs, with
30 epochs without validation loss improvement in every fit. All six axes
evaluated the same 462, 459 and 462 conditions in folds 0, 1 and 2.
The fold RBF bandwidths were 0.430414, 0.425762 and 0.439392. Every metric
and every disease score was defined in all three folds.

| Delta metric | HVG-DEG 3k | HVG-DEG 4k | HVG-DEG 5k | Denoise-DEG 3k | Denoise-DEG 4k | Denoise-DEG 5k |
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

| Signed disease DEG F1 | HVG-DEG 3k | HVG-DEG 4k | HVG-DEG 5k | Denoise-DEG 3k | Denoise-DEG 4k | Denoise-DEG 5k |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nHF | 0.1603 | 0.1639 | 0.1905 | 0.1197 | 0.1063 | 0.0969 |
| ICM | 0.1083 | 0.1258 | 0.1734 | 0.1211 | 0.1057 | 0.0961 |
| genetic cardiomyopathy | 0.2136 | 0.2469 | 0.2193 | 0.1208 | 0.1057 | 0.0950 |
| TOX | 0.3026 | 0.2587 | 0.2162 | 0.1302 | 0.1148 | 0.1036 |

See the [nine-metric figure](figure/delta-metrics.png) and the
[disease F1 figure](figure/disease-f1.png). The
[metric CSV](figure/data/delta-metrics.csv) and
[disease CSV](figure/data/disease-f1.csv) retain the full numerical
precision and eligibility counts. The materialized report retains the
per-fold and per-condition values and the full CM direction reference.

## Historical single-dose comparison

The [archive](historical/04_model_view/README.md) retains the accepted
single-dose comparison, receipts, payloads and original rendering data.
Use `python stage2/04_tsg_axis/execution/verify_historical.py --full` to
verify the sealed assets and replay all 54 fits' fixed predictions, the
nine aggregate artifacts and the three CSV-only figures. Replay writes to
scratch. The archive retains the original run plans, checkpoints and receipts.

## Replay the preceding axis-wide-control comparison

[Spec 5](https://github.com/Moke110/arc/issues/310) used an axis-wide
control for validation and test inputs and predicted Delta. It selected
the checkpoint by treated-expression MSE with 50-epoch patience. The
[original formal receipt](log/316-formal.json) records that comparison.
Its figures and plot CSVs are available in Git revision `85bc014`.
Its checkpoints, report and audit retain `output/formal`, `output/report`
and `log/316-runs`.

Use the Git checkout matching the registration when replaying archived
checkpoints. Prepare the axis runtime worktree under the checkout's root `tmp/`.

```bash
git worktree add --detach "$PWD/tmp/axis-runtime" 986f8ea
```

Use the same Python and CUDA environment as the recorded fit. The Step 04
dispatcher verifies the exact Git revision. For the 18 preceding formal fits,
pass the materialized matrix, Step 03 folder, frozen KPGT and native model source.

```bash
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
python stage2/04_tsg_axis/execution/replay_pinned.py axes \
  --runtime "$PWD/tmp/axis-runtime" --matrix stage2/04_tsg_axis/output/formal \
  --step03 stage2/03_data_train --representations stage2/util/output/representations \
  --source pert_model/TranSiGen/source
```

Replay the predecessor's scores and CSV-only figures with the runtime that
reads the separately materialized `log/316-runs` audit.

```bash
git worktree add --detach "$PWD/tmp/axis-report-runtime" 85bc014
STEP=stage2/04_tsg_axis
python "$PWD/tmp/axis-report-runtime/$STEP/execution/report_formal.py" verify \
  --matrix "$STEP/output/formal" --matrix-audit "$STEP/log/316-runs" \
  --report "$STEP/output/report" --figures "$PWD/tmp/axis-report-runtime/$STEP/figure" \
  --step03 stage2/03_data_train --raw stage2/01_data_raw/output \
  --logs "$PWD/tmp/axis-predecessor-report-verify"
```

## Replay downstream sealed runtimes

```bash
git worktree add --detach "$PWD/tmp/dose-width-runtime" e4587a8
```

For a recorded Step 07 or 08 checkpoint, use the same dispatcher with
`dose` or `width`, `--runtime "$PWD/tmp/dose-width-runtime"`, `--run` for
the sealed checkpoint folder, `--inputs` for the recorded tenv KPGT input
freeze and `--source` for the native TranSiGen source. Add `--scratch`
under `tmp/` for width replay. Step 07 checks exact validation and test
predictions. Step 08 checks exact validation predictions and held-out
test-score CSVs. The [compatibility receipt](log/310-runtime-replay.json)
records those checks against the accepted generations.
