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
the nine Step 04 scores on the common held-out rows.
