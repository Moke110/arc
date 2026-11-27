# Step 06.2 TranSiGen prediction architectures

Read the [architecture study routing](../AGENTS.md), the [study report](README.md)
and [ADR 0096](../../../docs/adr/0096-transigen-shape-round2.md). For runner or
scoring changes, use the [shared TranSiGen APIs](../../util/README.md#transigen-training-and-evaluation).

## Register and fit

Use the frozen tenv cold-drug folds, ordered Denoise-DEG 5k axis, measured
same-plate DMSO and seed 42. Materialize the [Step 03.7 training release](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json),
the frozen tenv KPGT representation release under `cpst/util/output/representations`
and the [read-only TranSiGen source](../../../pert_model/TranSiGen/manifest.sha256).
Preserve the 32 repeated training-fold standardized log-dose features and the
[Step 05.2 weight-1.75 protocol](../../05_tsg_loss/05.2_mse_l1/README.md)
except for the registered architecture and objective changes. Both direct
Δ readouts add the signed prediction to the measured control. In
`output-delta-nox2`, preserve the measured treated encoder and all three KL
terms, omit the treated reconstruction MSE and double the Δ MSE in training
and checkpoint selection. Keep the Step 06.2 FiLM head's identity
initialization and unnormalized modulation.

Use [`run_architecture.py`](execution/run_architecture.py) for individual
coordinates. Use [`matrix_architecture.py`](execution/matrix_architecture.py)
to register the entire roster before fitting. Run one CUDA fit at a time,
with prediction and independent replay in separate processes. Run from the
materialized repository root in `arc` and define `capped_step04` from the
[capped runtime](../../04_tsg_axis/AGENTS.md#prepare-the-capped-runtime).
For scoring, materialize the [Step 03.1 HVG-DEG 5k release](../../03_data_train/03.1_tenv_hvg-5k/output/manifest.json)
for fold bandwidth calibration. Follow the
[execution environment](../../../docs/agents/execution-environment.md).
Record WSL RAM, available Windows host RAM and GPU use every ten seconds.
Qualify runtime changes with the [bounded smoke](#reproduce-the-bounded-cuda-smoke).

## Reproduce the bounded CUDA smoke

Cover all eight arms on fold 0 with the registered limits. Keep the
coordinator sampler used in the formal route running while the matrix runs.
Require exact fit replay, independent scoring and byte-identical PNGs from
the two CSV-only renders. [Accepted smoke](README.md#reproduce-the-bounded-cuda-smoke)

```bash
set -e
step="$PWD/cpst/06_tsg_shape/06.2_round2"
mkdir -p "$PWD/tmp"
work=$(mktemp -d "$PWD/tmp/tsg-prediction-architecture-smoke.XXXXXX")
mkdir -p "$work/scratch" "$work/matplotlib"
export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"

capped_step04 python "$step/execution/matrix_architecture.py" --smoke \
  --output "$work/matrix" --logs "$work/audit"
for phase in score verify; do
  capped_step04 python "$step/execution/evaluate_architecture.py" "$phase" \
    --matrix "$work/matrix" --output "$work/report" --scratch "$work/scratch"
done
for output in figures figures-replay; do
  python "$step/execution/render_architecture.py" --smoke \
    --data "$work/report" --output "$work/$output"
done
diff -r "$work/figures" "$work/figures-replay"
```

## Reproduce the formal matrix

Reserve 8 GiB of scratch. Keep the candidate and scratch on one filesystem.
Reuse the three sealed Step 06.1 originals from
`06.1_round1/output/formal/prior-347/matrix`, with their report beside the
matrix. Pass `--prior-matrix` to `report_architecture.py` if that release
resides elsewhere. Require exact validation predictions, validation
objective and held-out predictions for the 21 new fits. Require parity
for the originals' fold scores and eligibility counts and run independent
scoring of the combined report. Pair each fold with its `original` score.
[Step 06.1 provenance](../06.1_round1/log/350-formal.json) ·
[Sealed predecessor](../06.1_round1/README.md#verify-the-accepted-347-source)

```bash
set -e
step="$PWD/cpst/06_tsg_shape/06.2_round2"
mkdir -p "$PWD/tmp"
work=$(mktemp -d "$PWD/tmp/tsg-prediction-architecture-formal.XXXXXX")
mkdir -p "$work/scratch" "$work/matplotlib"
export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
samples="$work/coordinator-resources.jsonl"

python "$step/../06.1_round1/execution/sample_resources.py" --output "$samples" &
sampler_pid=$!
stop_sampler() {
  kill -TERM "$sampler_pid" 2>/dev/null || true
  wait "$sampler_pid"
}
trap 'stop_sampler || true' EXIT
while [ ! -s "$samples" ]; do
  kill -0 "$sampler_pid" 2>/dev/null || exit 1
  sleep 1
done

capped_step04 python "$step/execution/matrix_architecture.py" \
  --output "$work/matrix" --logs "$work/matrix-audit"
# Keep sampling through the first reading after the matrix ends.
sleep 10
stop_sampler
trap - EXIT

capped_step04 python "$step/execution/report_architecture.py" build \
  --matrix "$work/matrix" --matrix-audit "$work/matrix-audit" \
  --report "$work/report" --publication "$work/publication" \
  --activity "$work/report-audit" --scratch "$work/scratch" \
  --coordinator-samples "$samples" --code-revision "$(git rev-parse HEAD)"
```

### Verify the installed #353 publication

Use the pinned `87564cda35c7887de2cf6c79d1b2f0221fc368e8` checkout, its
`stage2/` source layout and its `capped_step04` runtime. The published
verifier binds the original run code and historical plot files. Run the
following from that checkout with its accepted inputs materialized. Set
`release` to the accepted Round 2 release if its ignored payload is in
another checkout. Keep installed release files read-only and use a fresh
activity directory and scratch. [Pinned #353 recipe](https://github.com/Moke110/carps/blob/87564cda35c7887de2cf6c79d1b2f0221fc368e8/stage2/06_tsg_shape/06.2_round2/README.md#reproduce-the-formal-matrix) ·
[Accepted formal evidence](README.md#reproduce-the-formal-matrix)

```bash
set -e
step="$PWD/stage2/06_tsg_shape/06.2_round2"
release="$step/output/formal"
mkdir -p "$PWD/tmp"
work=$(mktemp -d "$PWD/tmp/tsg-prediction-architecture-verify.XXXXXX")
mkdir -p "$work/scratch" "$work/matplotlib"
export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
capped_step04 python "$step/execution/report_architecture.py" verify \
  --matrix "$release/matrix" --matrix-audit "$release/matrix-audit" \
  --report "$release/report" --publication "$step" \
  --activity "$work/report-audit" --scratch "$work/scratch"
```

## Replay the published figures

Run from the current repository root in the renderer environment recorded
in the receipt. Render into scratch and require byte-identical PNGs.
[Current figure evidence](README.md#replay-the-published-figures) ·
[Figure receipt](log/figure-render.json)

```bash
set -e
step="$PWD/cpst/06_tsg_shape/06.2_round2"
mkdir -p "$PWD/tmp"
work=$(mktemp -d "$PWD/tmp/tsg-architecture-figures.XXXXXX")
MPLCONFIGDIR="$work/matplotlib" python "$step/figure/script/render_comparison.py" \
  --data "$step/figure/data" --output "$work/figures"
for figure in learning-curve.png metrics.png metrics-paired.png; do
  cmp "$step/figure/$figure" "$work/figures/$figure"
done
```
