# TranSiGen prediction-axis comparison

[Spec 5](https://github.com/Moke110/arc/issues/310) compares three HVG-DEG
and three Denoise-DEG axes on the same tenv plate–parent–dose conditions.
Each fit starts with random weights at seed 42 and uses the original
TranSiGen internal widths, frozen 2,304-dimensional KPGT and the registered
latent numeric dose transform. The selected prediction axis sets the gene
interface width. Fold membership comes from the frozen cold-drug folds in
the [Step 03 tenv training release](../03_data_train/README.md).

Training uses the matched plate control and measured treated expression.
Validation and test use the overall reference for their axis. The best
checkpoint minimizes validation treated-expression MSE. At the fixed
prediction seam, `prediction.npy` stores absolute treated expression,
`predicted_delta.npy` is predicted expression minus the overall reference
and `measured_delta.npy` is measured treatment minus its matched plate
control. All three matrices follow `conditions.csv` and `genes.csv`.

## Run the formal matrix

The matrix driver checks the complete ordered condition and fold rosters
across all six releases. It registers all 18 coordinates before training,
then runs each fit, prediction and exact replay in separate processes.
For the formal dataset, use the full cold-drug folds and seed 42.

```bash
STEP=stage2/04_tsg_axis
WORK=$PWD/tmp/tsg-axis-formal
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
python "$STEP/execution/matrix.py" --output "$WORK/matrix"
python "$STEP/execution/report_formal.py" build \
  --matrix "$WORK/matrix" --report "$WORK/report" --figures "$WORK/figures" \
  --logs "$WORK/logs"
```

Use a single capped CUDA job for fitting. Run report generation and replay
after training completes. The accepted [formal receipt](log/316-formal.json)
pins all 18 run manifests, the report and its two CSV-reproducible PNGs.
The materialized paths are `output/formal` for the checkpoints, registrations,
predictions and per-process resource logs, and `output/report` for the
per-condition metrics, fold and eligibility tables, bandwidths and shared
disease reference. The two summary CSVs live next to the figure source.

```bash
python "$STEP/execution/report_formal.py" verify \
  --matrix "$STEP/output/formal" --report "$STEP/output/report" --figures "$STEP/figure" \
  --logs "$WORK/verify"
```

## Register and replay one fit

```bash
STEP=stage2/04_tsg_axis
RUN=$PWD/tmp/tsg-axis/hvg-deg-5k-fold-0
python "$STEP/execution/run.py" register --axis hvg-deg-5k --fold 0 \
  --training stage2/03_data_train/03.1_tenv_hvg-deg-5k/output --device cuda --output "$RUN"
python "$STEP/execution/run.py" train --output "$RUN"
python "$STEP/execution/run.py" predict --output "$RUN"
python "$STEP/execution/run.py" replay --output "$RUN"
```

Pass `--representations` and `--source` when frozen KPGT or immutable native
TranSiGen is materialized in a different checkout. Keep the same CPU thread
caps for each phase. Registration binds dataset and gene order, parent
and condition partitions, feature and native code digests and the runtime.
Prediction exactly replays validation before producing held-out predictions.
A fresh replay verifies all scores' source arrays and row identities against
the recorded values and publishes the run manifest. Repeating replay checks
the sealed run again without changing it. Pass the corresponding path flags
when an input has moved. Its registered digests still apply.

Pass `--smoke` only for a bounded two-epoch fit with 64 training, 24
validation and 24 test conditions. Its manifest marks it as smoke and it is
excluded from formal results. The [acceptance receipt](log/314-smoke-and-migration.json)
records the verified 5k fit and archived prediction replay.

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
retrospective results on published fixed axes built with tenv evidence.
The [evaluation acceptance receipt](log/315-metrics-smoke.json) records the
bounded evaluation and figure replay.

## Formal results

All 18 fits passed exact validation and held-out prediction replay. All six
axes evaluated the same 462, 459 and 462 conditions in folds 0, 1 and 2.
The fold RBF bandwidths were 0.430414, 0.425762 and 0.439392. Every metric
and every disease score was defined in all three folds.

| Delta metric | HVG-DEG 3k | HVG-DEG 4k | HVG-DEG 5k | Denoise-DEG 3k | Denoise-DEG 4k | Denoise-DEG 5k |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| JAC_top200 | 0.0602 | 0.0532 | 0.0526 | 0.0575 | 0.0555 | 0.0598 |
| JAC_pos100 | 0.0622 | 0.0402 | 0.0425 | 0.0492 | 0.0486 | 0.0481 |
| JAC_neg100 | 0.0631 | 0.0696 | 0.0664 | 0.0672 | 0.0646 | 0.0744 |
| Δ PCC | 0.0464 | 0.0399 | 0.0496 | 0.0410 | 0.0561 | 0.0571 |
| Δ Spearman | 0.0591 | 0.0627 | 0.0711 | 0.0453 | 0.0559 | 0.0577 |
| Δ RMSE | 0.4238 | 0.3865 | 0.3518 | 0.3819 | 0.3498 | 0.3143 |
| Δ E-distance | 0.0682 | 0.0915 | 0.0785 | 0.0575 | 0.0530 | 0.0500 |
| Δ MMD | 0.2882 | 0.3152 | 0.2753 | 0.2491 | 0.2280 | 0.2009 |
| Δ R² | -0.6011 | -0.7141 | -0.6488 | -0.7249 | -0.6680 | -0.6456 |

| Signed disease DEG F1 | HVG-DEG 3k | HVG-DEG 4k | HVG-DEG 5k | Denoise-DEG 3k | Denoise-DEG 4k | Denoise-DEG 5k |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| nHF | 0.1415 | 0.1472 | 0.1723 | 0.1056 | 0.0933 | 0.0859 |
| ICM | 0.0946 | 0.1118 | 0.1542 | 0.1049 | 0.0917 | 0.0841 |
| genetic cardiomyopathy | 0.1889 | 0.2210 | 0.1975 | 0.1063 | 0.0927 | 0.0841 |
| TOX | 0.2637 | 0.2283 | 0.1913 | 0.1128 | 0.0994 | 0.0901 |

See the [nine-metric figure](figure/delta-metrics.png) and the
[disease F1 figure](figure/disease-f1.png). The [metric CSV](figure/data/delta-metrics.csv)
and [disease CSV](figure/data/disease-f1.csv) retain the full numerical
precision and eligibility counts. The materialized report retains the
per-fold and per-condition values and the full CM direction reference.

## Historical single-dose comparison

The [archive](historical/04_model_view/README.md) retains the accepted
single-dose comparison, receipts, payloads and original rendering data.
Run `python stage2/04_tsg_axis/execution/verify_historical.py` against a
materialized archive to verify the sealed assets and replay an accepted
fixed prediction.
