# TranSiGen loss studies

## Selected objective

The selected TranSiGen loss is
`L = 0.1 * K + MSE_X1 + MSE_X2 + MSE_Δ + 1.75 * L1_Δ`.
`K` is the sum of the three native KL terms. Use measured same-plate DMSO
as model input `X1` and as the common reference for predicted and observed Δ.
The results assume a measured matched control is available.

This is the `alpha = 1.75` coordinate in Step 05.2, with unit weights on
both expression reconstruction MSEs and the Δ MSE. On the three frozen tenv
cold-drug folds on Denoise-DEG 5k, mean Δ PCC increased from 0.406610 to
0.431920 and mean Δ RMSE decreased from 0.259575 to 0.249215 compared with
the zero-weight reference. These are the best mean Δ PCC and Δ RMSE among
the published Step 05 coordinates. Weight 1.5 gave marginally lower
Δ E-distance and Δ MMD within Step 05.2.
[Step 05.2 equal-fold scores](05.2_mse_l1/figure/data/delta-metrics.csv) ·
[Stage 2 trade-offs](../README.md#active-experiments-and-releases)

## Loss studies

[Step 05.1](05.1_loss_metrics/README.md) compares MSE, cosine, Pearson correlation
and L1 Delta losses with coefficients 0.75, 1 and 1.25 on the denoise-5k axis.

[Step 05.2](05.2_mse_l1/README.md) varies the weight of an additive Delta
L1 term while retaining the native TranSiGen objective.

[Step 05.3](05.3_mse-l1_mse-x/README.md) varies the shared expression
reconstruction MSE weight at fixed Delta L1 weight 1.75.

[Step 05.4](05.4_mse-l1_pcc/README.md) and
[Step 05.5](05.5_mse-l1_cos/README.md) vary the weights of bounded
Delta PCC and cosine similarity losses, respectively, using the unit
reconstruction weight and fixed Delta L1 weight from Step 05.3.

## Verify the current figures

The [figure render receipt](log/343-figure-render.json) binds all 15 current
PNGs to the 15 plotted CSVs from the accepted study reports. The following command
compares two fresh CSV-only renders with the published PNGs and checks every
input hash against its original formal receipt. It uses the `arc` rendering
environment and requires neither model checkpoints nor prediction matrices.

```bash
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/tsg-loss-figure-check.XXXXXX")
python stage2/05_tsg_loss/figure/script/verify_figures.py --work "$WORK/replay"
```

Bars show the group mean, dots show the individual fold/seed runs, and error
bars show one sample standard deviation in each direction. Signed differences
pair each run with its matching reference before aggregation. Up and down
arrows identify the preferred direction of each metric. Learning curves show
one model per panel and one training and validation line per run.

## Reproduce the formal Delta similarity sweeps

Run from the materialized checkout with the `arc` environment and define
the capped launcher. [Step 04 capped runtime](../04_tsg_axis/README.md#prepare-the-capped-runtime)

The shared driver registers all 36 coordinates before training, runs one CUDA
fit at a time, and resumes from validated phase receipts at the same paths.
It samples WSL RAM and GPU usage every ten seconds and records peak RSS for
registration, training, prediction and replay.

Allow at least 6 GiB of free scratch space. Set `TSG_SIMILARITY_WORK_ROOT`
to an existing directory under this checkout's `tmp/` when using a backed
volume. Keep scratch, reports and the candidate publication on that volume
so the atomic publication rename stays on the same filesystem.

```bash
set -e
LOSS="$PWD/stage2/05_tsg_loss"
SHARED="$LOSS/similarity"
mkdir -p "$PWD/tmp"
WORK_ROOT="${TSG_SIMILARITY_WORK_ROOT:-$PWD/tmp}"
WORK=$(mktemp -d "$WORK_ROOT/tsg-similarity-formal.XXXXXX")
mkdir -p "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib"
REVISION=$(git rev-parse HEAD)
capped_step04 python -u "$SHARED/matrix_similarity.py" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
capped_step04 python "$SHARED/report_similarity.py" build \
  --matrix "$WORK/matrix" --matrix-audit "$WORK/matrix-audit" \
  --reports "$WORK/reports" --publication "$WORK/publication" \
  --activity "$WORK/report-audit" --scratch "$WORK/scratch" \
  --code-revision "$REVISION"
```

The driver accepts `--training`, `--representations` and `--source` to
locate pinned inputs in another checkout. The report accepts those options
plus `--bandwidth-training` and `--previous-step` for the sealed Step 05.3
release. The report prepares both studies before moving the verified
candidate into place. Repeat its command with `verify` in place of `build`
and omit `--code-revision` to check the publication again.
