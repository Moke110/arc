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

## Formal results

The comparison covers 45 coordinates across frozen cold-drug folds 0, 1
and 2 and training seeds 17, 42 and 2026. All 27 new fits exactly replayed
the validation objective, validation predictions and held-out predictions
in separate processes. The 18 references retain their accepted Step
06.1, 06.2 and 06.3 lineage. The evaluation used fixed validation and
test sampling seeds 22617 and 100042.

| Arm | Δ PCC mean ± SD | Δ RMSE mean ± SD |
| --- | ---: | ---: |
| `XX-XX` | 0.426340 ± 0.014345 | 0.250777 ± 0.002121 |
| `XX-XXΔ` | 0.427947 ± 0.019165 | 0.241382 ± 0.002937 |
| `X-XΔ` | 0.347833 ± 0.102598 | 0.248649 ± 0.007588 |
| `XΔ-XΔ` | 0.428830 ± 0.018574 | 0.240752 ± 0.002598 |
| `XΔ-XΔΔ` | 0.429835 ± 0.017268 | 0.240819 ± 0.003052 |

Each mean gives equal weight to the nine fold and seed coordinates. SD
is the sample standard deviation across coordinates. The independent
measured Δ decoder had the highest mean Δ PCC, and the shared measured
Δ decoder had the lowest mean Δ RMSE.

| Paired contrast | Δ PCC | Δ RMSE | Δ Spearman | Δ E-distance |
| --- | ---: | ---: | ---: | ---: |
| `XΔ-XΔΔ − XX-XXΔ` | +0.001888 ± 0.007564 | −0.000563 ± 0.001480 | +0.001995 ± 0.007217 | +0.001090 ± 0.001945 |
| `XΔ-XΔ − XΔ-XΔΔ` | −0.001005 ± 0.010334 | −0.000067 ± 0.000836 | −0.000546 ± 0.010250 | −0.000210 ± 0.002162 |

The prespecified comparisons subtract scores within each fold and seed
before averaging. Foldwise means give equal weight to the three training
repeats. For both contrasts, the SD across the nine pairs exceeds the
absolute mean difference on every focus metric. The held-out folds are
the three data splits, with three repeated training fits per split.

- [All nine metric summaries](figure/data/delta-metrics.csv) and
  [405 held-out values with eligibility counts](figure/data/fold-scores.csv)
- [72 paired values](figure/data/primary-contrasts.csv),
  [24 foldwise summaries](figure/data/primary-contrast-by-fold.csv) and
  [eight overall summaries](figure/data/primary-contrast-summary.csv)
- [Fit ledger](log/runs.csv), [epoch objectives](figure/data/epoch-objectives.csv)
  and [frozen fold bandwidths](log/bandwidths.csv)
- [Metric figures](figure/delta-metrics.png), [differences from the original](figure/delta-metrics-paired-difference.png),
  [paired contrasts](figure/primary-contrasts.png) and [learning curves](figure/learning-curves.png)
- [Formal acceptance receipt](log/359-formal.json), [matrix resources](log/359-resources.jsonl),
  [coordinator resources](log/359-coordinator-resources.jsonl) and
  [figure replay receipt](log/figure-render.json)

## Reproduce the formal study

Use the materialized checkout, `arc` environment and `capped_round4`
function defined above. Create a fresh work directory and verify the
accepted reference fits before starting any new fits. The coordinator
records GPU, WSL and Windows host availability every ten seconds through
the matrix.

```bash
set -e
STEP="$PWD/stage2/06_tsg_shape/06.4_round4"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-measured-encoder-formal.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"

capped_round4 python "$STEP/execution/prior_measured.py" > "$WORK/source-lineage.json"
python "$STEP/../06.1_round1/execution/sample_resources.py" \
  --output "$WORK/coordinator-resources.jsonl" &
SAMPLER_PID=$!
stop_sampler() {
  kill -TERM "$SAMPLER_PID" 2>/dev/null || true
  wait "$SAMPLER_PID"
}
trap 'stop_sampler || true' EXIT
until [ -s "$WORK/coordinator-resources.jsonl" ]; do
  kill -0 "$SAMPLER_PID" || exit 1
  sleep 1
done
capped_round4 python "$STEP/execution/matrix_measured.py" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
sleep 10
stop_sampler
trap - EXIT
capped_round4 python "$STEP/execution/publish_measured.py" build \
  --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
  --report "$WORK/report" --publication "$WORK/publication" \
  --activity "$WORK/publication-audit" --scratch "$WORK/scratch" \
  --coordinator-samples "$WORK/coordinator-resources.jsonl" --code-revision "$(git rev-parse HEAD)"
```

Use the same publisher with `verify` to check the installed release.
Provide the canonical `output/formal/{matrix,matrix-audit,report}`
paths, the step directory as `--publication` and fresh work directories
for `--activity` and `--scratch`. The verifier independently regenerates
the scientific scores, audits paired tables and replays the plots twice.
The figures also replay independently of scientific files through the
[Stage 2 figure verifier](../../README.md#verify-the-current-figures).
