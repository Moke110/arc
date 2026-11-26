# Step 06.2 TranSiGen prediction architectures

[ADR 0096](../../../docs/adr/0096-transigen-shape-round2.md) registers
eight prediction architectures on the Denoise-DEG 5k axis. Each fit uses
measured same-plate DMSO, frozen KPGT, and 32 repeats of log-dose
standardized on its training fold. The original expression encoder and
decoder widths are 1,200 and 800 with 100-dimensional latents.

| Architecture | Registered change |
| --- | --- |
| `original` | Step 06.1 original and Step 05.2 objective. |
| `output-delta` | Independent signed Δ decoder of the predicted treated latent. |
| `output-delta-nox2` | Use the treated decoder for signed Δ and omit the treated reconstruction. |
| `z-interact_MLP` | 532→145→200 prediction MLP with 106,485 parameters. |
| `z-interact_FiLM` | 432→136→200 conditioner and 100→200 projection with 106,488 parameters and identity modulation at initialization. |
| `shared-z` | Tied encoder trunks and separate Gaussian projections. |
| `z-res` | Sampled control latent enters the predicted Gaussian mean as a residual. |
| `shared-z-res` | Tied encoder trunks and a prediction-mean residual. |

Both direct Δ arms construct predicted treated expression as `X1 + Δ_hat`.
The `output-delta-nox2` training and validation objective is
`0.1 K + MSE_X1 + 2 MSE_Δ + 1.75 L1_Δ`. All three native KL terms remain
in `K`. The other arms use
`0.1 K + MSE_X1 + MSE_X2 + MSE_Δ + 1.75 L1_Δ`.
The linear residual changes the mean parameterization and is absorbable into
the projection weights.

## Reproduce the bounded CUDA smoke

The [acceptance receipt](log/352-smoke.json) records eight fold-0 runs with
seed 42, two epochs, and 64 train, 24 validation, and 24 test conditions in
the registered order. Train, validation and test use the matched controls.
Every fit has exact validation and test prediction replay. The score report
and all three figures are checked in separate processes.

Run from the materialized checkout in the `arc` environment. Use the
[Step 04 capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime).
Keep the coordinator resource sampler running while the matrix runs.

```bash
set -e
STEP="$PWD/stage2/06_tsg_shape/06.2_round2"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-prediction-architecture-smoke.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib" "$WORK/plot-csv"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
capped_step04 python "$STEP/execution/matrix_architecture.py" --smoke \
  --output "$WORK/matrix" --logs "$WORK/audit"
for phase in score verify; do
  capped_step04 python "$STEP/execution/evaluate_architecture.py" "$phase" \
    --matrix "$WORK/matrix" --output "$WORK/report" --scratch "$WORK/scratch"
done
cp "$WORK/report/"{delta-metrics,fold-scores,epoch-objectives}.csv "$WORK/plot-csv"
for output in figures figures-replay; do
  python "$STEP/figure/script/render_architecture.py" --smoke \
    --data "$WORK/plot-csv" --output "$WORK/$output"
done
diff -r "$WORK/figures" "$WORK/figures-replay"
```
