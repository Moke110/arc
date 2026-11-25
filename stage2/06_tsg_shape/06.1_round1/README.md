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

## Formal results

All 24 fits completed independent exact validation and test replay. The
three fresh `widen_z-100` fits reproduced the accepted Step 05.2 weight-1.75
training and validation traces, selected checkpoint tensors, validation and
test predictions, Delta arrays, and all 27 fold scores and eligibility
counts. The selected epochs were 129, 116, and 85 for folds 0, 1, and 2.

[Mean held-out scores](figure/delta-metrics.png) ·
[Mean paired differences](figure/delta-metrics-paired-difference.png) ·
[Learning curves](figure/learning-curves.png) ·
[216 fold scores](figure/data/fold-scores.csv) ·
[72 equal-fold aggregates](figure/data/delta-metrics.csv) ·
[Epoch objectives](figure/data/epoch-objectives.csv)

Bars show equal-fold means for both the held-out scores and signed
differences from the fresh control. Dots show the individual replicates.
Error bars show sample SD within each architecture. Panel arrows mark
the preferred metric direction.

Each learning-curve panel shows one architecture with train and validation
objectives for each fold and seed through the last epoch observed for
that fit. Objectives are displayed on log y-scales. The wider encoder has
its own y-scale to show its full objective range.

| Coordinate | JAC_top200 | JAC_pos100 | JAC_neg100 | Δ PCC | Δ Spearman | Δ RMSE | Δ E-distance | Δ MMD | Δ R² |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `widen_z-100` | 0.152470 | 0.129491 | 0.183801 | 0.431920 | 0.413273 | 0.249215 | 0.025826 | 0.154295 | -0.009312 |
| `widen_z-110` | 0.149212 | 0.126256 | 0.180178 | 0.425745 | 0.407999 | 0.251052 | 0.027748 | 0.161386 | -0.043429 |
| `widen_z-120` | 0.149755 | 0.126661 | 0.180951 | 0.425208 | 0.408778 | 0.251897 | 0.027044 | 0.159553 | -0.044712 |
| `widen_z-130` | 0.151713 | 0.128085 | 0.183928 | 0.424276 | 0.406629 | 0.251211 | 0.026343 | 0.157022 | -0.037053 |
| `widen_enc` | 0.151442 | 0.126972 | 0.183873 | 0.426133 | 0.410575 | 0.250727 | 0.029191 | 0.165253 | -0.032227 |
| `deepen_enc` | 0.148972 | 0.125434 | 0.180330 | 0.426294 | 0.409516 | 0.249230 | 0.024126 | 0.147206 | -0.022935 |
| `deepen_dec` | 0.152206 | 0.139344 | 0.171715 | 0.431106 | 0.417729 | 0.249401 | 0.021348 | 0.140121 | 0.005436 |
| `deepen-both` | 0.149561 | 0.134907 | 0.171445 | 0.428210 | 0.414693 | 0.249063 | 0.020687 | 0.134288 | -0.000714 |

Increasing the latent width gave lower mean Δ PCC and higher mean Δ RMSE
at all three tested widths. Deepening the decoder raised JAC_pos100 from
0.129491 to 0.139344 and lowered JAC_neg100 from 0.183801 to 0.171715.
Deepening both sides lowered Δ E-distance from 0.025826 to 0.020687 and
Δ MMD from 0.154295 to 0.134288. These results describe the three frozen
tenv folds at seed 42.

The widest encoder used 71,640,600 trainable parameters and reached about
1.41 GiB of CUDA tensor allocations. WSL retained at least 17.29 GiB of
available RAM during the 23-minute matrix. The [formal receipt](log/347-formal.json)
binds the matrix and exact replay evidence. The [run ledger](log/runs.csv),
[KL diagnostics](log/kl-diagnostics.csv), [matrix resource samples](log/347-resources.jsonl),
and [coordinator samples](log/347-coordinator-resources.jsonl) give the
per-run and workstation evidence. The accepted local artifacts are under
`output/formal/{matrix,matrix-audit,report}`.

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
  python "$STEP/figure/script/render_shape.py" --smoke --data "$WORK/report" --output "$WORK/$target"
done
diff -r "$WORK/figures" "$WORK/figures-replay"
```

## Reproduce the formal matrix

Use the materialized checkout and the same capped runtime as the smoke.
Reserve at least 8 GiB of scratch space. When necessary, use a backed
volume reachable under the checkout’s `tmp/` directory. Keep the candidate
publication and report scratch on that volume so the publication rename
stays on one filesystem. Keep a coordinator GPU and RAM sample stream on
the same monotonic clock at ten-second intervals throughout the matrix.

```bash
set -e
STEP="$PWD/stage2/06_tsg_shape/06.1_round1"
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

## Replay the published figures

The three canonical figures use the replicate, aggregate and epoch CSVs. Use
the versions and digests in the [figure receipt](log/figure-render.json).
The formal replay command verifies the scientific receipt and compares the
current CSV-only renders against the figure receipt.

```bash
STEP="$PWD/stage2/06_tsg_shape/06.1_round1"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-shape-figures.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
python "$STEP/figure/script/render_shape.py" \
  --data "$STEP/figure/data" --output "$WORK/figures"
for figure in delta-metrics.png delta-metrics-paired-difference.png learning-curves.png; do
  cmp "$STEP/figure/$figure" "$WORK/figures/$figure"
done
```
