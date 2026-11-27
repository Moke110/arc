# Step 06.4 TranSiGen measured encoder targets

The numbered commands use the [shared TranSiGen training and evaluation APIs](../../util/README.md#transigen-training-and-evaluation).

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
evidence. To reproduce the figure hashes in that receipt, use renderer
revision `f387f2c3`. The current renderer uses a distinct hue for each
arm in the metric figures under [ADR 0093](../../../docs/adr/0093-stage2-shared-plotting-and-in-place-figures.md).

Run from the materialized checkout with the `arc` environment active.
Start the [coordinator sampler](../06.1_round1/execution/sample_resources.py)
before the matrix and continue it for at least one sample after completion.
Check for 12 GiB of available WSL memory before each capped command. Use
the commands below with an unused work directory.

```bash
set -e
STEP="$PWD/cpst/06_tsg_shape/06.4_round4"
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
stop_sampler() {
  kill -TERM "$SAMPLER_PID" 2>/dev/null || true
  wait "$SAMPLER_PID"
}
trap 'stop_sampler || true' EXIT
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
  -s "$STEP/../../util/training/tests/06.4_round4" -p 'test_measured_model.py'
stop_sampler
trap - EXIT
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
- [Metric figures](figure/metrics.png), [differences from the original](figure/metrics-paired.png),
  [learning curves](figure/learning-curve.png)
- [Original paired contrast figure](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.4_round4/figure/primary-contrasts.png)
- [Formal acceptance receipt](log/359-formal.json), [matrix resources](log/359-resources.jsonl),
  [coordinator resources](log/359-coordinator-resources.jsonl) and
  [figure replay receipt](log/figure-render.json)

## Reproduce the formal study

The [original release reproduction](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.4_round4/README.md#reproduce-the-formal-study) pins the renderer and publication verifier for the historical results. Use that revision for the formal release verification below. The current comparison figures have their own CSV-only replay.

Use the materialized checkout, `arc` environment and `capped_round4`
function defined above. Create a fresh work directory and verify the
accepted reference fits before starting any new fits. The coordinator
records GPU, WSL and Windows host availability every ten seconds through
the matrix.

```bash
set -e
STEP="$PWD/cpst/06_tsg_shape/06.4_round4"
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

## Replay a sealed Step 06.4 fit

The [formal receipt](log/359-formal.json) pins the matrix code at
`f387f2c370eb87cd27c02502aa98f165851491a7`. Its registrations bind code by
repository-relative path as well as digest. Use a detached worktree with
the original `stage2/` layout. Supply the canonical locations of the same
materialized inputs. Run the following from the live checkout with
`capped_round4` defined above. The example copies one sealed shared decoder
fit into scratch and replays that copy.

```bash
BASE="$PWD"
mkdir -p "$BASE/tmp"
WORK=$(mktemp -d "$BASE/tmp/tsg-measured-pinned.XXXXXX")
STAGE="$BASE/cpst"
STEP="$STAGE/06_tsg_shape/06.4_round4"
SEALED_RUN="$STEP/output/formal/matrix/runs/XΔ-XΔ/fold-0/seed-17"
git worktree add --detach "$WORK/code" f387f2c370eb87cd27c02502aa98f165851491a7
mkdir -p "$WORK/scratch" "$WORK/matplotlib" "$WORK/run"
cp --reflink=auto --dereference -- "$SEALED_RUN"/* "$WORK/run/"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib" PYTHONDONTWRITEBYTECODE=1
capped_round4 python "$WORK/code/stage2/06_tsg_shape/06.4_round4/execution/run_measured.py" replay \
  --output "$WORK/run" \
  --training "$STAGE/03_data_train/03.7_tenv_denoise-5k/output" \
  --representations "$STAGE/util/output/representations" \
  --source "$BASE/pert_model/TranSiGen/source"
git worktree remove "$WORK/code"
```

## Replay the published figures

The three current figures use the [epoch traces](figure/data/learning-curve.csv),
[held-out scores](figure/data/metrics.csv) and [paired differences](figure/data/metrics-paired.csv).
Differences subtract `XX-XX` within each fold and seed. The
[replay receipt](log/figure-render.json) records the accepted fit evidence,
plotted CSV and PNG digests, and renderer versions.

```bash
STEP="$PWD/cpst/06_tsg_shape/06.4_round4"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-architecture-figures.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
python "$STEP/figure/script/render_comparison.py" \
  --data "$STEP/figure/data" --output "$WORK/figures"
for figure in learning-curve.png metrics.png metrics-paired.png; do
  cmp "$STEP/figure/$figure" "$WORK/figures/$figure"
done
```
