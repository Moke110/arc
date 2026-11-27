# TranSiGen prediction-axis comparison

Read the [accepted comparison](README.md) and follow the
[CarPS-Transcriptome workflow](../AGENTS.md). For runner or scoring changes,
use the [shared TranSiGen APIs](../util/README.md#transigen-training-and-evaluation).
Study-specific input binding, metrics and disease reference code live in
[`tsg_axis`](output/api/tsg_axis).

## Prepare the capped runtime

Activate `arc` and work from the repository root with the releases materialized.
Fit on the six 3k, 4k and 5k HVG-DEG and Denoise-DEG plate-pooled tenv releases
with the [shared five-fold tenv evaluation release](../util/README.md#shared-tenv-evaluation). Use the Step 02 ordered axis tables, frozen
2,304-dimensional tenv KPGT vectors at `cpst/util/output/representations/`, and
native source at `pert_model/TranSiGen/source/`. Scoring also needs the Step 01
raw tenv corpus and the Step 02 HVG-DEG disease reference. Follow the owner
routes for input preparation. [Training releases](../03_data_train/README.md) ·
[Gene axes](../02_gene_axis/README.md) · [KPGT reader](../util/README.md#shared-input-and-release-apis) ·
[TranSiGen provenance](../../pert_model/TranSiGen/author_log.json) · [Raw corpus](../01_data_raw/README.md)

Run full fits serially and monitor available RAM and GPU usage every 10 seconds
while fitting. Use the [execution environment](../../docs/agents/execution-environment.md)
for resource limits and bounded reads. Keep candidates under `tmp/`.

```bash
step="$PWD/cpst/04_tsg_axis"
mkdir -p "$PWD/tmp"

capped_step04() {
  if ! LC_ALL=C free -g | awk '
    { print }
    $1 == "Mem:" { available = $7 }
    END { exit available < 12 }
  '; then
    printf '%s\n' 'Step 04 requires at least 12 GiB of available RAM.' >&2
    return 1
  fi
  systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- \
    env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 "$@"
}
```

## Run the formal matrix

[`matrix.py`](execution/matrix.py) checks the ordered condition and fold rosters
across the six releases and registers all 150 axis, fold and seed coordinates before training,
and runs train, predict and replay in separate processes. Require exact validation
objective, validation predictions and held-out predictions for each fit.
[`report_formal.py`](execution/report_formal.py) checks the matrix and process
receipts, scores and verifies the report, and checks two independent PNG renders.

```bash
work=$(mktemp -d "$PWD/tmp/tsg-axis-matrix.XXXXXX")
matrix=(python "$step/execution/matrix.py"
  --output "$work/matrix" --logs "$work/matrix-audit/logs")
capped_step04 "${matrix[@]}" --register-only &&
capped_step04 "${matrix[@]}" &&
capped_step04 python "$step/execution/report_formal.py" build \
  --matrix "$work/matrix" --matrix-audit "$work/matrix-audit" \
  --report "$work/report" --figures "$work/figures" --logs "$work/report-logs"
```

Run the same matrix command after an interruption with the same inputs and
output directories. The matrix authenticates the existing registrations
against the planned axis, fold, seed and partition before fitting any
coordinate. It continues from verified phase receipts. All 150 are new
registrations because the pinned three-fold runs use different memberships.

The report binds its per-fold full-axis bandwidth to the shared fold-release
digest and the actual training rows. Its `bandwidths.csv` is the Step 09
calibration lineage input. Install accepted successors under
`output/formal-tenv-v1/` and `output/report-tenv-v1/` while preserving the
historical releases described below. Keep the process logs beside the matrix
in `output/formal-tenv-v1-audit/logs/` and keep the coordinator samples with
that audit. Verify the installed result before replacing the tracked figures.

```bash
work=$(mktemp -d "$PWD/tmp/tsg-axis-verification.XXXXXX")
capped_step04 python "$step/execution/report_formal.py" verify \
  --matrix "$step/output/formal-tenv-v1" \
  --matrix-audit "$step/output/formal-tenv-v1-audit" \
  --report "$step/output/report-tenv-v1" --figures "$step/figure" \
  --logs "$work/report-logs"
```

Each registered fit records the full runtime and source hashes. Use the
source commit and environment in the formal receipt to resume or replay a
fit when the active training code has changed.

The study uses matched same-plate DMSO as the model input in every partition.
At the fixed-prediction boundary, `prediction.npy` is absolute treated expression.
Subtract the same matched control from prediction and measured treatment to form
both Δ arrays in registered condition and gene order. Keep held-out treated
measurements for scoring. Select the checkpoint by the lowest per-condition
native validation objective (MSE + KLD, β = 0.1), with 30 epochs of patience
and a 1,000-epoch maximum. [Matched-control protocol](../../docs/adr/0086-tenv-axis-matched-plate-control.md)

## Register and replay one fit

This example runs the bounded HVG-DEG 5k smoke through exact replay and report
verification. `--smoke` bounds the fit to two epochs and 64 train, 24 validation
and 24 test conditions. For a full single-coordinate fit, omit `--smoke` at
registration and stop after `replay`. Formal scoring requires the full matrix.
[Bounded gate evidence](log/399-tenv-v1-formal.json)

```bash
work=$(mktemp -d "$PWD/tmp/tsg-axis-smoke.XXXXXX")
run="$work/run"
report="$work/report"
capped_step04 python "$step/execution/run.py" register --axis hvg-5k --fold 4 --seed 2027 \
  --training cpst/03_data_train/03.1_tenv_hvg-5k/output --device cuda --smoke --output "$run" &&
capped_step04 python "$step/execution/run.py" train --output "$run" &&
capped_step04 python "$step/execution/run.py" predict --output "$run" &&
capped_step04 python "$step/execution/run.py" replay --output "$run" &&
capped_step04 python "$step/execution/evaluate.py" score --smoke --run "$run" --output "$report" &&
capped_step04 python "$step/execution/evaluate.py" verify --smoke --run "$run" --output "$report"
```

For relocated pinned inputs, pass `--training`, `--representations` or `--source`
to the fit phases, and `--step03` or `--raw` to the evaluator. Keep the registered
runtime and input digests. For scoring changes, run
`python -m unittest discover -s cpst/util/eval/tests/04_tsg_axis`.
[Runner](execution/run.py) · [Evaluator](execution/evaluate.py) ·
[Fixed-prediction checks](../util/eval/fixed_prediction.py)

## Verify the sealed report

The installed `output/formal-matched-control/` runs and
`output/report-matched-control/` retain their historical axis IDs and sealed
digests. `log/319-runs/` holds the process logs. Verify with the evaluator at
`3f4cbb65` and alias the six Step 03 paths to the axis names it expects. The
pinned worktree carries its disease-reference tables. Copy the compact report to
scratch for the verification build. [Formal receipt](log/319-formal.json) ·
[Sealed evaluator and axis mapping](log/322-canonical-axis-report.json)

```bash
(
  set -e
  work=$(mktemp -d "$PWD/tmp/tsg-axis-verify.XXXXXX")
  git worktree add --detach "$work/runtime" 3f4cbb65
  trap 'git worktree remove "$work/runtime"' EXIT
  mkdir "$work/legacy-step03"
  for release in "$PWD"/cpst/03_data_train/03.*_tenv_*-[345]k; do
    name=${release##*/}
    ln -s "$release" "$work/legacy-step03/${name%-*}-deg-${name##*-}"
  done
  runs=()
  for run in "$step"/output/formal-matched-control/runs/*/fold-*; do
    runs+=(--run "$run")
  done
  cp -a "$step/output/report-matched-control" "$work/report"
  capped_step04 python "$work/runtime/stage2/04_tsg_axis/execution/evaluate.py" verify \
    "${runs[@]}" --output "$work/report" \
    --step03 "$work/legacy-step03" --raw "$PWD/cpst/01_data_raw/output"
)
```

For canonical IDs in the six axis-bearing tables, run
[`export_canonical.py`](execution/export_canonical.py). It checks the sealed
digests and preserves every non-axis field.

```bash
work=$(mktemp -d "$PWD/tmp/tsg-axis-export.XXXXXX")
python "$step/execution/export_canonical.py" --output "$work/report"
```

## Rebuild and verify the comparison figures

Use `report_formal.py` to reproduce the successor's plotted tables and PNGs
from its installed matrix. The paired figure uses `denoise-5k` as the
within-fold and seed reference. The
[rendering receipt](log/figure-render.json) binds the active plots and tables.
For the pinned three-fold matrix and its report, the historical
[`build_comparison.py`](execution/build_comparison.py) still checks
JAC-top200 parity against the sealed condition and fold scores.

For figure-only verification, render from the published CSVs and compare the
three PNGs byte for byte. [CSV-only renderer](figure/script/render.py)

```bash
(
  set -e
  work=$(mktemp -d "$PWD/tmp/tsg-axis-render.XXXXXX")
  MPLCONFIGDIR="$work/mpl" python "$step/figure/script/render.py" \
    --data "$step/figure/data" --output "$work/replayed"
  for name in learning-curve.png metrics.png metrics-paired.png; do
    cmp "$step/figure/$name" "$work/replayed/$name"
  done
)
```

When checking independence from the retired archive, run
[`verify_active_inputs.py`](execution/verify_active_inputs.py) through
`capped_step04` with Steps 01–03 materialized, including Tahoe, in a checkout
without `cpst/history`. Check for `PASS` and `history_available: false` in the
JSON result. [Archive-isolated smoke receipt](log/323-archive-isolated-smoke.json)
