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

## Register and replay a fit

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

## Historical single-dose comparison

The [archive](historical/04_model_view/README.md) retains the accepted
single-dose comparison, receipts, payloads and original rendering data.
Run `python stage2/04_tsg_axis/execution/verify_historical.py` against a
materialized archive to verify the sealed assets and replay an accepted
fixed prediction.
