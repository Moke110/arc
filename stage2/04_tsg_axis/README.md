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

## Historical single-dose comparison

The [archive](historical/04_model_view/README.md) retains the accepted
single-dose comparison, receipts, payloads and original rendering data.
Run `python stage2/04_tsg_axis/execution/verify_historical.py` against a
materialized archive to verify the sealed assets and replay an accepted
fixed prediction.
