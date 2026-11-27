# Step 06.3 TranSiGen architecture combinations

Read the [architecture study routing](../AGENTS.md), the [study report](README.md)
and [ADR 0097](../../../docs/adr/0097-transigen-shape-round3.md). For runner or
scoring changes, use the [shared TranSiGen APIs](../../util/README.md#transigen-training-and-evaluation).

## Register and fit

Materialize the [Step 03.7 training release](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json),
the frozen tenv KPGT release under `cpst/util/output/representations` and
the [read-only TranSiGen source](../../../pert_model/TranSiGen/manifest.sha256).
Keep the [Step 06.2 protocol](../../../docs/adr/0096-transigen-shape-round2.md)
for the three cold-drug folds, ordered Denoise-DEG 5k axis, measured same-plate
DMSO and 32 repeats of training-fold standardized log dose. Use training seeds
17, 42, 2026 and 43 with validation and test sampling seeds 22617 and 100042.

In arms 3–5, preserve FiLM's identity initialization and apply
`LayerNorm(100, elementwise_affine=False)` after modulation and before the
Gaussian projection. For direct Δ readouts, add the signed prediction to the
measured control. Arms 2 and 5 retain the measured treated encoder and all
three KL terms and replace treated reconstruction MSE with a second Δ MSE term.

Use [`run_composition.py`](execution/run_composition.py) for individual fits and
[`matrix_composition.py`](execution/matrix_composition.py) to register the
new fits before training. Run one CUDA fit at a time with separate training,
prediction and independent replay processes. Apply the
[execution environment](../../../docs/agents/execution-environment.md) and
[capped runtime](../../04_tsg_axis/README.md#prepare-the-capped-runtime)
from the materialized repository root in the `arc` environment.
Record available WSL and Windows host RAM and GPU use every ten seconds.
Materialize the [Step 03.1 HVG-DEG 5k release](../../03_data_train/03.1_tenv_hvg-5k/output/manifest.json)
for fold bandwidth calibration before scoring.

## Reproduce the bounded CUDA smoke

Use all six arms on folds 0 and 1 and seeds 17 and 42 for two epochs.
The matrix caps each fold at the first 64 training, 24 validation and 24
held-out conditions. Require exact replay for all 24 fits, an objective audit,
independent scoring verification and two matching CSV-only renders.
[Initial and revised smoke evidence](README.md#first-attempt-and-normalization)

```bash
(
  set -e
  STEP="$PWD/cpst/06_tsg_shape/06.3_round3"
  mkdir -p "$PWD/tmp"
  WORK=$(mktemp -d "$PWD/tmp/tsg-architecture-composition-smoke.XXXXXX")
  mkdir -p "$WORK/scratch" "$WORK/matplotlib" "$WORK/plot-csv"
  export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
  SAMPLES="$WORK/coordinator-resources.jsonl"

  python "$STEP/../06.1_round1/execution/sample_resources.py" --output "$SAMPLES" &
  sampler_pid=$!
  trap 'kill -TERM "$sampler_pid" 2>/dev/null || true; wait "$sampler_pid" || true' EXIT
  while [ ! -s "$SAMPLES" ]; do
    kill -0 "$sampler_pid" 2>/dev/null || exit 1
    sleep 1
  done
  capped_step04 python "$STEP/execution/matrix_composition.py" --smoke \
    --output "$WORK/matrix" --logs "$WORK/matrix-audit"
  kill -TERM "$sampler_pid"
  wait "$sampler_pid"
  trap - EXIT

  python "$STEP/execution/audit_smoke_objectives.py" --matrix "$WORK/matrix"
  python "$STEP/execution/prior_composition.py" > "$WORK/source-lineage.json"
  for phase in score verify; do
    capped_step04 python "$STEP/execution/evaluate_composition.py" "$phase" --smoke \
      --matrix "$WORK/matrix" --output "$WORK/report" --scratch "$WORK/scratch"
  done
  cp "$WORK/report/"{delta-metrics,fold-scores,epoch-objectives,primary-contrasts,primary-contrast-summary}.csv \
    "$WORK/plot-csv"
  for output in figures figures-replay; do
    python "$STEP/execution/render_composition.py" --smoke \
      --data "$WORK/plot-csv" --output "$WORK/$output"
  done
  diff -r "$WORK/figures" "$WORK/figures-replay"
)
```

To verify the original sealed smoke, run `report_composition.py verify --smoke`
at `7ce9e6259c6554072a7193b2830ce9fd03c54248` with its original matrix, matrix
audit, report and publication. Its matrix ran at
`29c7e1680e79b32d40de18950260ff57ad5615bd`.
[Original smoke procedure](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.3_round3/README.md#reproduce-the-bounded-cuda-smoke)

## Reproduce the revised formal run

Complete the [bounded smoke](#reproduce-the-bounded-cuda-smoke) before starting
the formal matrix. Keep the nine sealed seed-42 fits from
[Rounds 1 and 2](../06.2_round2/log/353-formal.json) and the 27 reusable fits
in `output/formal/prior-unnormalized` from the [first Round 3 attempt](log/356-blocked.json)
read-only. Record at least one coordinator sample before the matrix starts
and one after it finishes. Keep the candidate and report scratch on the same
filesystem. If the reusable fits have moved, pass their verified locations
with `--round-one-matrix`, `--round-two-matrix` and `--prior-unnormalized`.

```bash
(
  set -e
  STEP="$PWD/cpst/06_tsg_shape/06.3_round3"
  mkdir -p "$PWD/tmp"
  WORK=$(mktemp -d "$PWD/tmp/tsg-architecture-composition-formal.XXXXXX")
  mkdir -p "$WORK/scratch" "$WORK/matplotlib"
  export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
  REVISION=$(git rev-parse HEAD)
  SAMPLES="$WORK/coordinator-resources.jsonl"

  python "$STEP/../06.1_round1/execution/sample_resources.py" --output "$SAMPLES" &
  sampler_pid=$!
  trap 'kill -TERM "$sampler_pid" 2>/dev/null || true; wait "$sampler_pid" || true' EXIT
  while [ ! -s "$SAMPLES" ]; do
    kill -0 "$sampler_pid" 2>/dev/null || exit 1
    sleep 1
  done
  capped_step04 python "$STEP/execution/matrix_composition.py" \
    --output "$WORK/matrix" --logs "$WORK/matrix-audit"
  sleep 10
  kill -TERM "$sampler_pid"
  wait "$sampler_pid"
  trap - EXIT

  capped_step04 python "$STEP/execution/report_composition.py" build \
    --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
    --report "$WORK/report" --publication "$WORK/publication" \
    --activity "$WORK/report-audit" --scratch "$WORK/scratch" \
    --coordinator-samples "$SAMPLES" --code-revision "$REVISION"
)
```

The matrix must seal exact validation objective, validation prediction and
held-out prediction replays for the 36 new fits. The report checks source
artifacts, lineage and ordered populations before reuse and independently
rescores all 72 coordinates. Pair each score with `0_original` within the
same fold and seed and retain the ADR 0097 component and interaction contrasts.
[Formal receipt](log/356-formal.json) · [Table audit](log/356-formal-table-audit.json)

### Verify the installed historical release

Use `28f5eea4593170466ff6ca4110def1ef3b125cce`, which retains the
`stage2/` paths bound by the sealed registrations and the accepted report
verifier. Match the recorded runtime in the installed registrations. The
verifier rescans the stored fits and regenerates the historical publication
for comparison. Copy the report to scratch so its publication lock stays
outside the installed release.
[Installed release inventory](../../util/training/verification.md) ·
[Original formal publication](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.3_round3/README.md#reproduce-the-revised-formal-run)

```bash
(
  set -e
  ROOT="$PWD"
  STAGE="$ROOT/cpst"
  SHAPE="$STAGE/06_tsg_shape"
  FORMAL="$SHAPE/06.3_round3/output/formal"
  mkdir -p "$ROOT/tmp"
  WORK=$(mktemp -d "$ROOT/tmp/tsg-architecture-composition-verify.XXXXXX")
  mkdir -p "$WORK/scratch" "$WORK/matplotlib"
  export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"

  git worktree add --detach "$WORK/runtime" 28f5eea4593170466ff6ca4110def1ef3b125cce
  HISTORICAL="$WORK/runtime/stage2/06_tsg_shape/06.3_round3"
  cp -a --reflink=auto "$FORMAL/report" "$WORK/report"
  capped_step04 python "$HISTORICAL/execution/report_composition.py" verify \
    --matrix "$FORMAL/matrix" --matrix-audit "$FORMAL/matrix-audit" \
    --report "$WORK/report" --publication "$HISTORICAL" \
    --activity "$WORK/report-audit" --scratch "$WORK/scratch" \
    --training "$STAGE/03_data_train/03.7_tenv_denoise-5k/output" \
    --bandwidth-training "$STAGE/03_data_train/03.1_tenv_hvg-5k/output" \
    --representations "$STAGE/util/output/representations" \
    --source "$ROOT/pert_model/TranSiGen/source" \
    --round-one-matrix "$SHAPE/06.1_round1/output/formal/prior-347/matrix" \
    --round-two-matrix "$SHAPE/06.2_round2/output/formal/matrix" \
    --prior-unnormalized "$FORMAL/prior-unnormalized"
)
```

## Replay the published figures

Render the three current PNGs from the published plot CSVs in the software
versions recorded in the figure receipt. [Figure inputs and receipt](README.md#replay-the-published-figures)

```bash
(
  set -e
  STEP="$PWD/cpst/06_tsg_shape/06.3_round3"
  mkdir -p "$PWD/tmp"
  WORK=$(mktemp -d "$PWD/tmp/tsg-architecture-composition-figures.XXXXXX")
  mkdir -p "$WORK/scratch" "$WORK/matplotlib"
  export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
  python "$STEP/figure/script/render_comparison.py" \
    --data "$STEP/figure/data" --output "$WORK/figures"
  for figure in learning-curve.png metrics.png metrics-paired.png; do
    cmp "$STEP/figure/$figure" "$WORK/figures/$figure"
  done
)
```
