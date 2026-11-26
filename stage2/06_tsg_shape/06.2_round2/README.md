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

## Formal results

On the three frozen tenv cold-drug folds, FiLM and the prediction-mean
residual each increased mean Δ PCC and reduced mean Δ RMSE from the
original. Both direct Δ readouts reduced Δ RMSE on every fold and had
lower mean Δ PCC. The comparison includes the readout changes in both arms, plus the
omitted treated reconstruction and doubled Δ MSE in `output-delta-nox2`.

| Architecture | Δ PCC | Δ RMSE | JAC_top200 | JAC_pos100 | JAC_neg100 | Δ Spearman | Δ E-distance | Δ MMD | Δ R² |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `original` | 0.431920 | 0.249215 | 0.152470 | 0.129491 | 0.183801 | 0.413273 | 0.025826 | 0.154295 | -0.009312 |
| `output-delta` | 0.423164 | 0.241952 | 0.160278 | 0.144541 | 0.183792 | 0.400691 | 0.026695 | 0.137706 | 0.087260 |
| `output-delta-nox2` | 0.415948 | 0.241974 | 0.158069 | 0.142596 | 0.181327 | 0.393649 | 0.028720 | 0.141267 | 0.085782 |
| `z-interact_MLP` | 0.423109 | 0.252523 | 0.150532 | 0.126175 | 0.182706 | 0.408493 | 0.029128 | 0.165812 | -0.061340 |
| `z-interact_FiLM` | 0.439568 | 0.247354 | 0.161090 | 0.137472 | 0.192909 | 0.421419 | 0.025326 | 0.149580 | -0.022400 |
| `shared-z` | 0.408751 | 0.254451 | 0.139689 | 0.118574 | 0.168399 | 0.394348 | 0.026461 | 0.156450 | -0.159086 |
| `z-res` | 0.436771 | 0.246442 | 0.154476 | 0.137375 | 0.179482 | 0.417866 | 0.020049 | 0.131439 | 0.011988 |
| `shared-z-res` | 0.424692 | 0.247712 | 0.148903 | 0.131022 | 0.174834 | 0.407311 | 0.020195 | 0.130342 | 0.000678 |

FiLM improved both primary metrics on two folds, and `z-res` reduced
Δ RMSE on all three. The similarly sized MLP and the shared encoder alone
had lower Δ PCC and higher Δ RMSE on every fold. Adding the residual to
the shared encoder improved both primary means relative to `shared-z`.
Compared with `original`, the combined arm had lower mean Δ PCC and lower
mean Δ RMSE. All fits used seed 42 and require measured matched DMSO at
prediction time.

[Held-out scores](figure/delta-metrics.png) ·
[Paired differences](figure/delta-metrics-paired-difference.png) ·
[Learning curves](figure/learning-curves.png) ·
[216 fold scores and eligibility counts](figure/data/fold-scores.csv) ·
[72 aggregates and signed differences](figure/data/delta-metrics.csv) ·
[Epoch objectives](figure/data/epoch-objectives.csv)

The score and paired-difference panels show equal-fold means, individual
fold dots and sample SD. Each paired difference subtracts the original
score on the same fold. The learning curves identify each arm's objective
and show the train and validation traces for all three folds. The no-X2
objective uses its own log-scale y-axis.

All 21 new fits passed exact independent validation prediction, validation
objective and held-out prediction replay. The three Step 06.1 originals
reproduced all 27 baseline scores and eligibility counts. The formal report
was rescored independently, and its three PNGs were reproduced twice
from the CSVs. The [formal receipt](log/353-formal.json) binds these checks,
the [fit ledger](log/runs.csv), and the [matrix](log/353-resources.jsonl)
and [coordinator](log/353-coordinator-resources.jsonl) resource records.
The matrix finished in 14 minutes 54 seconds. Coordinator samples reached
a maximum interval of 10.002 seconds, with at least 17.12 GiB of WSL RAM
and 20.21 GiB of Windows host RAM available.
The local matrix, phase audits and scored report are under
`output/formal/{matrix,matrix-audit,report}`.

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

## Reproduce the formal matrix

Run from the materialized checkout in the `arc` environment and use the
[Step 04 capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime).
Reserve 8 GiB of scratch. Keep the publication candidate and its scratch
on the same filesystem. The accepted Step 06.1 originals remain in
`06.1_round1/output/formal/prior-347/matrix`. Supply `--prior-matrix` to
`report_architecture.py` if that verified source has moved.

```bash
set -e
STEP="$PWD/stage2/06_tsg_shape/06.2_round2"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-prediction-architecture-formal.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
REVISION=$(git rev-parse HEAD)
SAMPLES="$WORK/coordinator-resources.jsonl"

python "$STEP/../06.1_round1/execution/sample_resources.py" --output "$SAMPLES" &
SAMPLER_PID=$!
stop_sampler() {
  kill -TERM "$SAMPLER_PID" 2>/dev/null || true
  wait "$SAMPLER_PID"
}
trap 'stop_sampler || true' EXIT
while [ ! -s "$SAMPLES" ]; do
  kill -0 "$SAMPLER_PID" 2>/dev/null || exit 1
  sleep 1
done

capped_step04 python "$STEP/execution/matrix_architecture.py" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
# Leave the sampler active for a reading after the matrix ends.
sleep 10
stop_sampler
trap - EXIT

capped_step04 python "$STEP/execution/report_architecture.py" build \
  --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
  --report "$WORK/report" --publication "$WORK/publication" \
  --activity "$WORK/report-audit" --scratch "$WORK/scratch" \
  --coordinator-samples "$SAMPLES" --code-revision "$REVISION"
```

The driver registers all 21 coordinates before training. It fits one
coordinate at a time and runs prediction and exact validation and test
replay in separate processes. The report validates the three sealed
originals and independently reproduces all scores and PNGs. To verify an
installed release, run `report_architecture.py verify` with the installed
`output/formal/{matrix,matrix-audit,report}` and the step directory as
`--publication`. Omit `--coordinator-samples` and `--code-revision`, which
are taken from the publication. Use fresh paths for `--activity` and
`--scratch` under the checkout's `tmp/`.

## Replay the published figures

The renderer takes only the plotted CSVs. The
[figure receipt](log/figure-render.json) records their identities and the
tool versions for byte-identical replay.

```bash
STEP="$PWD/stage2/06_tsg_shape/06.2_round2"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-prediction-architecture-figures.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
python "$STEP/figure/script/render_architecture.py" \
  --data "$STEP/figure/data" --output "$WORK/figures"
for figure in delta-metrics.png delta-metrics-paired-difference.png learning-curves.png; do
  cmp "$STEP/figure/$figure" "$WORK/figures/$figure"
done
```
