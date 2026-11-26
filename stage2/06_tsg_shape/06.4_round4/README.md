# Step 06.4 TranSiGen measured encoder targets

[ADR 0098](../../../docs/adr/0098-transigen-shape-round4.md) compares five
architectures on the frozen tenv Denoise-DEG 5k protocol with measured
same-plate DMSO, KPGT and a numeric dose standardized on the training fold.
The stochastic prediction head uses the control latent and drug–dose
features in each arm.

| Arm | Measured encoders | Prediction readout |
| --- | --- | --- |
| `XX-XX` | Control and treated expression | Treated-expression decoder |
| `XX-XXΔ` | Control and treated expression | Independent signed Δ decoder |
| `X-XΔ` | Control | Signed Δ decoder |
| `XΔ-XΔ` | Control and measured signed Δ | Signed Δ decoder shared with reconstruction |
| `XΔ-XΔΔ` | Control and measured signed Δ | Independent signed Δ decoder |

Each direct Δ arm forms treated prediction as measured control plus
decoded Δ. The control-only arm uses the control Gaussian KL and
reconstruction, double Δ MSE and 1.75 times Δ L1. Both measured Δ arms
align the predicted Gaussian to the Gaussian encoding of measured
`X2 − X1` and reconstruct that signed measurement. The trace columns
`treatment_mse` and `treatment_kl` denote the measured branch's
reconstruction and KL. They contain signed Δ terms in the measured Δ arms
and zeros in the control-only arm.

## Bounded CUDA gate

The gate trains and independently replays all five architectures on the
first two frozen folds and seeds 17 and 42, with 64 training, 24
validation and 24 held-out conditions per fold for two epochs. It uses
the same fit, prediction, replay and report entry points as the formal
study. The [gate receipt](log/358-gate.json) records the matrix, objective,
source provenance, independent rescoring, CSV figure renders and resource
evidence.

Run from the materialized checkout with the `arc` environment active.
Start the [coordinator sampler](../06.1_round1/execution/sample_resources.py)
before the matrix and continue it for at least one sample after completion.
Check for 12 GiB of available WSL memory before each capped command. Use
the commands below with an unused work directory.

```bash
STEP="$PWD/stage2/06_tsg_shape/06.4_round4"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-measured-encoder-smoke.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"

capped_round4() {
  LC_ALL=C free -g | awk '$1 == "Mem:" {available=$7} END {exit available<12}' || return
  systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
    env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 "$@"
}

python "$STEP/../06.1_round1/execution/sample_resources.py" \
  --output "$WORK/coordinator-resources.jsonl" &
SAMPLER_PID=$!
until [ -s "$WORK/coordinator-resources.jsonl" ]; do
  kill -0 "$SAMPLER_PID" || exit 1
  sleep 1
done
capped_round4 python "$STEP/execution/matrix_measured.py" --smoke \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
sleep 10
capped_round4 python "$STEP/execution/audit_measured.py" objectives --smoke --matrix "$WORK/matrix"
for phase in score verify; do
  capped_round4 python "$STEP/execution/evaluate_measured.py" "$phase" --smoke \
    --matrix "$WORK/matrix" --output "$WORK/report" --scratch "$WORK/scratch"
done
python "$STEP/execution/prior_measured.py"
python "$STEP/execution/audit_measured_scores.py" --smoke --report "$WORK/report"
python "$STEP/execution/audit_measured_figures.py" build --smoke \
  --data "$WORK/report" --output "$WORK/figures" --receipt "$WORK/figures.json" \
  --scratch "$WORK/scratch"
capped_round4 env TSG_MEASURED_TEST_DEVICE=cuda python -m unittest discover \
  -s "$STEP/execution/tests" -p 'test_measured_model.py'
kill -TERM "$SAMPLER_PID"
wait "$SAMPLER_PID"
python "$STEP/execution/audit_measured.py" resources --smoke \
  --matrix-audit "$WORK/matrix-audit" --coordinator "$WORK/coordinator-resources.jsonl" \
  --coverage-output "$WORK/coordinator-coverage.jsonl"
```

The formal roster uses all three folds and training seeds 17, 42 and 2026.
It retains 18 reference fits from the accepted Round 3 release and requires
27 new fits with exact independent validation-objective, validation-prediction
and held-out prediction replay.
