# Step 06.1 TranSiGen expression shape round one

Eight expression architectures share the ordered Denoise-DEG 5k gene axis,
frozen tenv cold-drug folds, seed 42, frozen KPGT vectors, standardized
numeric dose, and measured same-plate DMSO throughout. The objective is
`0.1 × K + MSE_X1 + MSE_X2 + MSE_Δ + 1.75 × L1_Δ`, where `K` sums the
three native KL terms. The earliest minimum validation objective determines
the checkpoint, with patience 30 and a cap of 1,000 epochs. Both expression
branches use the indicated widths and retain independent parameters. The
treated decoder serves both reconstruction and prediction.

| Coordinate | Encoder hidden widths | Decoder hidden widths | Expression latent dimension |
| --- | --- | --- | --- |
| `widen_z-100` | 1200 | 800 | 100 |
| `widen_z-110` | 1200 | 800 | 110 |
| `widen_z-120` | 1200 | 800 | 120 |
| `widen_z-130` | 1200 | 800 | 130 |
| `widen_enc` | 6000 | 800 | 100 |
| `deepen_enc` | 1200, 1200 | 800 | 100 |
| `deepen_dec` | 1200 | 800, 800 | 100 |
| `deepen-both` | 1200, 1200 | 800, 800 | 100 |

The fresh `widen_z-100` scores are the within-fold references. Equal-fold
means use defined fold scores. Mean signed differences use folds where both
scores are defined. The total KL diagnostic sums the three raw KL means per
condition. The per-latent-dimension diagnostic divides that sum by the
shared expression latent width. [Design](../../../docs/adr/0094-transigen-shape-round1.md)
 · [Step 05.2 protocol](../../05_tsg_loss/05.2_mse_l1/README.md)

## Reproduce the bounded smoke

The [CUDA smoke receipt](log/346-smoke.json) records all six short runs and
independent replays against the published training release. All eight shapes
were exercised by the model build test. The smoke uses the first 64 training,
24 validation, and 24 test conditions from fold 0 for two epochs.

Run from the materialized checkout in the `arc` environment. Define
`capped_step04` as shown in the [capped runtime setup](../../04_tsg_axis/README.md#prepare-the-capped-runtime).
Set `WORK` to a fresh directory below this checkout's `tmp/`, using a backed
volume if needed for the checkpoints.

```bash
set -e
STEP="$PWD/stage2/06_tsg_shape/06.1_round1"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-shape-smoke.XXXXXX")
mkdir -p "$WORK/tmp" "$WORK/matplotlib" "$WORK/runs"
export TMPDIR="$WORK/tmp" MPLCONFIGDIR="$WORK/matplotlib"
SHAPES=(widen_z-100 widen_z-130 widen_enc deepen_enc deepen_dec deepen-both)
RUNS=()
for shape in "${SHAPES[@]}"; do
  run="$WORK/runs/$shape"
  capped_step04 python "$STEP/execution/run_shape.py" register \
    --output "$run" --shape "$shape" --fold 0 --smoke
  RUNS+=(--run "$run")
done
for shape in "${SHAPES[@]}"; do
  for phase in train predict replay; do
    capped_step04 python "$STEP/execution/run_shape.py" "$phase" \
      --output "$WORK/runs/$shape"
  done
done
for phase in score verify; do
  capped_step04 python "$STEP/execution/evaluate_shape.py" "$phase" \
    --smoke "${RUNS[@]}" --output "$WORK/report" --scratch "$WORK/tmp"
done
for target in figures figures-replay; do
  python "$STEP/figure/script/render_shape.py" --data "$WORK/report" --output "$WORK/$target"
done
diff -r "$WORK/figures" "$WORK/figures-replay"
```
