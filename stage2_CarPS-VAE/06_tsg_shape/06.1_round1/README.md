# Step 06.1 TranSiGen expression shape round one

Ten expression architectures share the ordered Denoise-DEG 5k gene axis,
frozen tenv cold-drug folds, seed 42, frozen KPGT vectors, standardized
numeric dose, and measured same-plate DMSO throughout. The objective is
`0.1 × K + MSE_X1 + MSE_X2 + MSE_Δ + 1.75 × L1_Δ`, where `K` sums the
three native KL terms. The earliest minimum validation objective determines
the checkpoint, with patience 30 and a cap of 1,000 epochs. Both expression
branches use the indicated widths and retain independent parameters. The
treated decoder serves both reconstruction and prediction.

| Coordinate | Encoder hidden widths | Decoder hidden widths | Expression latent dimension |
| --- | --- | --- | --- |
| `original` | 1200 | 800 | 100 |
| `widen_z-110` | 1200 | 800 | 110 |
| `widen_z-120` | 1200 | 800 | 120 |
| `widen_z-130` | 1200 | 800 | 130 |
| `widen_enc-2k` | 2000 | 800 | 100 |
| `widen_enc-4k` | 4000 | 800 | 100 |
| `widen_enc-6k` | 6000 | 800 | 100 |
| `deepen_enc` | 1200, 1200 | 800 | 100 |
| `deepen_dec` | 1200 | 800, 800 | 100 |
| `deepen_both` | 1200, 1200 | 800, 800 | 100 |

The accepted `original` scores are the within-fold references. Equal-fold
means use defined fold scores. Mean signed differences use folds where both
scores are defined. The total KL diagnostic sums the three raw KL means per
condition. The per-latent-dimension diagnostic divides that sum by the
shared expression latent width. [Round-one design](../../../docs/adr/0094-transigen-shape-round1.md)
 · [Width extension](../../../docs/adr/0095-transigen-shape-encoder-width-extension.md)
 · [Step 05.2 protocol](../../05_tsg_loss/05.2_mse_l1/README.md)

## Formal results

The six new fits at widths 2,000 and 4,000 completed independent exact
validation and test replay across all three folds. The comparison reuses
the 24 accepted [#347](https://github.com/Moke110/carps/issues/347) fits. Their
three original-architecture fits reproduced the Step 05.2 weight-1.75
traces, selected checkpoints, validation and test predictions, Delta
arrays, and all 27 fold scores and eligibility counts. They selected
epochs 129, 116, and 85 for folds 0, 1, and 2.

[Mean held-out scores](figure/metrics.png) ·
[Mean paired differences](figure/metrics-paired.png) ·
[Learning curves](figure/learning-curve.png) ·
[270 fold scores](figure/data/fold-scores.csv) ·
[90 equal-fold aggregates](figure/data/delta-metrics.csv) ·
[Epoch objectives](figure/data/epoch-objectives.csv)

Bars show equal-fold means for both the held-out scores and signed
differences from `original` on the same fold. Dots show the individual replicates.
Error bars show sample SD within each architecture. Panel arrows mark
the preferred metric direction.

Each learning-curve panel shows one architecture with train and validation
objectives for each fold and seed through the last epoch observed for
that fit. Objectives are displayed on log y-scales. Each widened encoder
has its own y-scale to show the full range of its training spikes.

| Coordinate | JAC_top200 | JAC_pos100 | JAC_neg100 | Δ PCC | Δ Spearman | Δ RMSE | Δ E-distance | Δ MMD | Δ R² |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `original` | 0.152470 | 0.129491 | 0.183801 | 0.431920 | 0.413273 | 0.249215 | 0.025826 | 0.154295 | -0.009312 |
| `widen_z-110` | 0.149212 | 0.126256 | 0.180178 | 0.425745 | 0.407999 | 0.251052 | 0.027748 | 0.161386 | -0.043429 |
| `widen_z-120` | 0.149755 | 0.126661 | 0.180951 | 0.425208 | 0.408778 | 0.251897 | 0.027044 | 0.159553 | -0.044712 |
| `widen_z-130` | 0.151713 | 0.128085 | 0.183928 | 0.424276 | 0.406629 | 0.251211 | 0.026343 | 0.157022 | -0.037053 |
| `widen_enc-2k` | 0.152057 | 0.128353 | 0.183962 | 0.428968 | 0.412196 | 0.251051 | 0.029489 | 0.166656 | -0.055261 |
| `widen_enc-4k` | 0.153355 | 0.131377 | 0.183310 | 0.433229 | 0.416627 | 0.248576 | 0.027660 | 0.159373 | -0.016795 |
| `widen_enc-6k` | 0.151442 | 0.126972 | 0.183873 | 0.426133 | 0.410575 | 0.250727 | 0.029191 | 0.165253 | -0.032227 |
| `deepen_enc` | 0.148972 | 0.125434 | 0.180330 | 0.426294 | 0.409516 | 0.249230 | 0.024126 | 0.147206 | -0.022935 |
| `deepen_dec` | 0.152206 | 0.139344 | 0.171715 | 0.431106 | 0.417729 | 0.249401 | 0.021348 | 0.140121 | 0.005436 |
| `deepen_both` | 0.149561 | 0.134907 | 0.171445 | 0.428210 | 0.414693 | 0.249063 | 0.020687 | 0.134288 | -0.000714 |

At 4,000 encoder units, mean Δ PCC rose from 0.431920 to 0.433229 and
mean Δ RMSE fell from 0.249215 to 0.248576. At 2,000 and 6,000 units,
mean Δ PCC fell and mean Δ RMSE rose. Mean Δ E-distance and Δ MMD rose
at all three wider encoder settings.

Increasing the latent width gave lower mean Δ PCC and higher mean Δ RMSE
at all three tested widths. Deepening the decoder raised JAC_pos100 from
0.129491 to 0.139344 and lowered JAC_neg100 from 0.183801 to 0.171715.
Deepening both sides lowered Δ E-distance from 0.025826 to 0.020687 and
Δ MMD from 0.154295 to 0.134288. These results describe the three frozen
tenv folds at seed 42.

The new widths used 30,016,600 and 50,828,600 trainable parameters and
peaked at 0.64 and 1.03 GiB of CUDA tensor allocations, respectively.
WSL retained at least 17.55 GiB of available RAM in the six-minute matrix.
The [formal receipt](log/350-formal.json) binds both fit sources and their
exact replay evidence. The [repeat check](log/350-formal-repeat.json)
compares the first formal candidate with the six post-smoke fits and
records identical scores, selected epochs and learning curves. The [run ledger](log/runs.csv),
[KL diagnostics](log/kl-diagnostics.csv), [matrix resource samples](log/350-resources.jsonl),
and [coordinator samples](log/350-coordinator-resources.jsonl) give the
per-run and workstation evidence. The six-fit matrix, its audit and the
combined report are under `output/formal/{matrix,matrix-audit,report}`.

## Verify the accepted #347 source

The sealed matrix, audit and report remain under `output/formal/prior-347/`.
The [original receipt](log/347-formal.json) is pinned at the [#347 release revision](https://github.com/Moke110/carps/tree/9dee5fd90ea6d2bed775655b0205875c8812938b/stage2/06_tsg_shape/06.1_round1).
The three historical names `widen_z-100`, `widen_enc` and `deepen-both`
map to `original`, `widen_enc-6k` and `deepen_both` in the combined report.
The historical registrations, receipts and checkpoints keep their source
names and hashes. In the `arc` environment, use
`python stage2_CarPS-VAE/06_tsg_shape/06.1_round1/execution/report_shape.py verify-prior`
to verify them. Supply `--prior-matrix` if the source has moved, and keep its
`report` directory beside the matrix.

For CUDA replay of a historical run, use the training code pinned at
`1f166d2f0dee53444131eee1dd83481c8dedce4f` in a detached worktree. Pass the
original coordinate path to that revision's `execution/run_shape.py replay`.
Run the following from a materialized checkout with the capped runtime
available.

```bash
BASE="$PWD"
mkdir -p "$BASE/tmp"
WORK=$(mktemp -d "$BASE/tmp/tsg-shape-prior.XXXXXX")
STEP=stage2_CarPS-VAE/06_tsg_shape/06.1_round1
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib" PYTHONDONTWRITEBYTECODE=1
git worktree add --detach "$WORK/code" 1f166d2f0dee53444131eee1dd83481c8dedce4f
capped_step04 python "$WORK/code/stage2/06_tsg_shape/06.1_round1/execution/run_shape.py" replay \
  --output "$BASE/$STEP/output/formal/prior-347/matrix/runs/widen_enc/fold-0" \
  --training "$BASE/stage2_CarPS-VAE/03_data_train/03.7_tenv_denoise-5k/output" \
  --representations "$BASE/stage2_CarPS-VAE/util/output/representations" \
  --source "$BASE/pert_model/TranSiGen/source"
git worktree remove "$WORK/code"
```

## Reproduce the bounded smoke

The [CUDA smoke receipt](log/349-smoke.json) records the two added widths
and their independent replays against the published training release.
The [bounded report receipt](log/349-smoke-report.json) records the original
architecture's exact replay, all 27 scores for the three bounded fits, and
two byte-identical renders of both score figures. This check followed the
initial formal candidate and passed before the six fits in the published
matrix. The model build tests exercised all ten shapes. The smoke uses
the first 64 training, 24 validation, and 24 test conditions from fold 0
for two epochs. The commands below reproduce the complete bounded report
with the original architecture and both added widths.

Run from the materialized checkout in the `arc` environment. Define
`capped_step04` as shown in the [capped runtime setup](../../04_tsg_axis/README.md#prepare-the-capped-runtime).
Set `WORK` to a fresh directory below this checkout's `tmp/`, using a backed
volume if needed for the checkpoints.

```bash
set -e
STEP="$PWD/stage2_CarPS-VAE/06_tsg_shape/06.1_round1"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-shape-smoke.XXXXXX")
mkdir -p "$WORK/tmp" "$WORK/matplotlib" "$WORK/runs"
export TMPDIR="$WORK/tmp" MPLCONFIGDIR="$WORK/matplotlib"
SHAPES=(original widen_enc-2k widen_enc-4k)
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
  python "$STEP/execution/render_shape.py" --smoke --data "$WORK/report" --output "$WORK/$target"
done
diff -r "$WORK/figures" "$WORK/figures-replay"
```

## Reproduce the formal matrix

The [original release reproduction](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/06_tsg_shape/06.1_round1/README.md#reproduce-the-formal-matrix) pins the renderer and publication verifier for the historical results. Use that revision for the formal release verification below. The current comparison figures have their own CSV-only replay.

Use the materialized checkout and the same capped runtime as the smoke.
The driver registers, fits and replays only the six new coordinates. It
uses the sealed #347 source when building the combined report.
Reserve at least 8 GiB of scratch space. When necessary, use a backed
volume reachable under the checkout’s `tmp/` directory. Keep the candidate
publication and report scratch on that volume so the publication rename
stays on one filesystem. Keep a coordinator GPU and RAM sample stream on
the same monotonic clock at ten-second intervals throughout the matrix.

```bash
set -e
STEP="$PWD/stage2_CarPS-VAE/06_tsg_shape/06.1_round1"
mkdir -p "$PWD/tmp"
WORK_ROOT="${TSG_SHAPE_WORK_ROOT:-$PWD/tmp}"
WORK=$(mktemp -d "$WORK_ROOT/tsg-shape-formal.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
REVISION=$(git rev-parse HEAD)
COORDINATOR_SAMPLES="$WORK/coordinator-resources.jsonl"
python "$STEP/execution/sample_resources.py" --output "$COORDINATOR_SAMPLES" &
COORDINATOR_PID=$!
stop_coordinator() {
  kill -TERM "$COORDINATOR_PID" 2>/dev/null || true
  wait "$COORDINATOR_PID"
}
trap 'stop_coordinator || true' EXIT
while [ ! -s "$COORDINATOR_SAMPLES" ]; do
  kill -0 "$COORDINATOR_PID" 2>/dev/null || exit 1
  sleep 1
done

capped_step04 python "$STEP/execution/matrix_shape.py" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
stop_coordinator
trap - EXIT

capped_step04 python "$STEP/execution/report_shape.py" build \
  --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
  --prior-matrix "$STEP/output/formal/prior-347/matrix" \
  --report "$WORK/report" --publication "$WORK/publication" \
  --activity "$WORK/report-audit" --scratch "$WORK/scratch" \
  --coordinator-samples "$COORDINATOR_SAMPLES" --code-revision "$REVISION"
```

The coordinator sampler writes GPU 0 and WSL RAM every ten seconds on
the matrix's monotonic clock. It also samples Windows host RAM when
PowerShell is available. Start it in the same WSL distro before the
matrix and give it access to `nvidia-smi`. If a matrix run resumes, keep
the sampler running across attempts in that WSL boot to preserve coverage.
The matrix keeps its own resource samples and per-phase peak RSS. To
independently replay a completed coordinate again, invoke
`run_shape.py replay --output` on the coordinate directory.

`report_shape.py verify` repeats formal scoring, control parity and CSV
render verification. Use the same arguments as `build`, with the sealed
`output/formal/{matrix,matrix-audit,report}` and the step directory as
`--publication`. Omit `--coordinator-samples` and `--code-revision`.
The matrix accepts alternate locations for the verified training,
representations and upstream source. The publication command accepts
alternate verified training, bandwidth-training and Step 05.2 paths.
Keep the candidate separate from the accepted release through verification.

## Replay the published figures

The three current figures use the [epoch traces](figure/data/learning-curve.csv),
[held-out scores](figure/data/metrics.csv) and [paired differences](figure/data/metrics-paired.csv).
Differences subtract `original` within each fold and seed. The
[replay receipt](log/figure-render.json) records the accepted fit evidence,
plotted CSV and PNG digests, and renderer versions.

```bash
STEP="$PWD/stage2_CarPS-VAE/06_tsg_shape/06.1_round1"
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
