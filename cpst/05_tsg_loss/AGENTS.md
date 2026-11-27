# Step 05 TranSiGen loss studies

Read the [Step 05 report](README.md) for objective selection and follow the
[CarPS-Transcriptome workflow](../AGENTS.md). For runner or scoring changes,
use the [shared TranSiGen APIs](../util/README.md#transigen-training-and-evaluation).

## Study entry points

- [05.1](05.1_loss_metrics/AGENTS.md) compares Delta loss forms and coefficients.
- [05.2](05.2_mse_l1/AGENTS.md) sweeps the additive Delta L1 weight.
- [05.3](05.3_mse-l1_mse-x/AGENTS.md) sweeps the shared expression reconstruction MSE weight.
- [05.4](05.4_mse-l1_pcc/AGENTS.md) and [05.5](05.5_mse-l1_cos/AGENTS.md) share the bounded Delta similarity runner and own their PCC and cosine studies.

Replay accepted fits from scratch copies using the revision and paths pinned
by each study's registration.
[05.1 replay](05.1_loss_metrics/AGENTS.md#replay-a-sealed-fit) ·
[05.2 replay](05.2_mse_l1/AGENTS.md#replay-a-sealed-fit) ·
[05.3 replay](05.3_mse-l1_mse-x/AGENTS.md#replay-a-sealed-fit) ·
[PCC replay](05.4_mse-l1_pcc/AGENTS.md#replay-a-sealed-fit) ·
[Cosine replay](05.5_mse-l1_cos/AGENTS.md#replay-a-sealed-fit)

Keep study receipts under the numbered study's `log/` directory. Both
similarity smoke receipts pin the same resource capture in Step 05.4.
[Shared smoke resource capture](05.4_mse-l1_pcc/log/338-coordinator-resources.jsonl)

## Verify the current figures

Activate `arc` and run from the repository root. The verifier checks the
acceptance and render receipts, then makes two independent renders of all
15 canonical PNGs from the tracked plot CSVs. Require both renders to
match the published bytes. [CSV-only verifier](figure/script/verify_figures.py)

```bash
(
  set -e
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-loss-figure-check.XXXXXX")
  python cpst/05_tsg_loss/figure/script/verify_figures.py --work "$work/replay"
)
```

## Rebuild the current comparisons from accepted fits

Activate `arc` at the materialized repository root and define
`capped_step04`. Use the frozen tenv Denoise-DEG 5k training release and
the five accepted Step 05 fit and report releases. Keep those inputs
read-only. [Capped runtime](../04_tsg_axis/AGENTS.md#prepare-the-capped-runtime) ·
[Training manifest](../03_data_train/03.7_tenv_denoise-5k/output/manifest.json) ·
[Comparison builder](execution/build_comparison.py)

The builder checks accepted manifest and fit identities and the complete
train and validation traces. It verifies that fixed predictions and their
matched DMSO controls reproduce the sealed Delta arrays, and checks sealed
JAC-top200 parity before adding wJAC-2000. It writes plot CSVs, verified
PNGs and receipts to a fresh candidate. Pass `--evidence-root` for a
different materialized checkout, or `--step` with a numbered study
directory name to rebuild just that study.

```bash
(
  set -e
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-loss-rebuild.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
  capped_step04 python cpst/05_tsg_loss/execution/build_comparison.py \
    --output "$work/candidate"
)
```

The prior figure set and original publication commands are pinned in Git.
[Historical Step 05 publication](https://github.com/Moke110/carps/tree/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/05_tsg_loss)

## Reproduce the original formal Delta similarity sweeps

Use `arc` in the materialized repository root and define `capped_step04`.
Fit on tenv Denoise-DEG 5k with the frozen cold-drug folds 0–2, seed 42,
frozen tenv KPGT vectors and the pinned TranSiGen source. Use measured
same-plate DMSO as input and as the common reference for predicted and
observed Delta. Scoring uses the frozen HVG-DEG 5k bandwidths and checks
each fresh zero-weight fit against the accepted Step 05.3 unit-weight fit.
[Capped runtime](../04_tsg_axis/AGENTS.md#prepare-the-capped-runtime) ·
[Frozen folds](../03_data_train/output/frozen-folds/tenv/manifest.json) ·
[KPGT reader and pins](../util/README.md#shared-input-and-release-apis) ·
[TranSiGen source manifest](../../pert_model/TranSiGen/manifest.sha256) ·
[Step 05.3 receipt](05.3_mse-l1_mse-x/log/336-formal.json) ·
[Similarity protocol](../../docs/adr/0092-transigen-bounded-delta-similarity-sweeps.md)

The matrix registers all 36 coordinates before training and runs one CUDA
fit at a time. Resume with the same matrix and audit paths to reuse
validated phase receipts. The driver samples WSL RAM and GPU use every
ten seconds and records phase peak RSS. Allow at least 6 GiB of free
scratch space. For a backed volume, set `TSG_SIMILARITY_WORK_ROOT` to an
existing directory or backed symlink under the checkout's `tmp/`. Keep
scratch, reports and candidate publication on that volume for the atomic
publication rename. The report checks exact replay and both studies'
frozen references, checks two independent renders of each original
figure set, and verifies the candidate after publication.
[Matrix driver](similarity/matrix_similarity.py) ·
[Report builder and verifier](similarity/report_similarity.py)

```bash
(
  set -e
  shared="$PWD/cpst/05_tsg_loss/similarity"
  mkdir -p "$PWD/tmp"
  work_root="${TSG_SIMILARITY_WORK_ROOT:-$PWD/tmp}"
  work=$(mktemp -d "$work_root/tsg-similarity-formal.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
  revision=$(git rev-parse HEAD)
  capped_step04 python -u "$shared/matrix_similarity.py" \
    --output "$work/matrix" --logs "$work/matrix-audit"
  report_args=(
    --matrix "$work/matrix" --matrix-audit "$work/matrix-audit"
    --reports "$work/reports" --publication "$work/publication"
    --activity "$work/report-audit" --scratch "$work/scratch"
  )
  capped_step04 python "$shared/report_similarity.py" build "${report_args[@]}" \
    --code-revision "$revision"
  capped_step04 python "$shared/report_similarity.py" verify "${report_args[@]}"
)
```

Pass `--training`, `--representations` and `--source` to the matrix and
report commands when pinned inputs are materialized elsewhere. The report
also accepts `--bandwidth-training` for the HVG-DEG 5k release and
`--previous-step` for the sealed Step 05.3 study root. Require both study
receipts to pass.
[PCC acceptance](05.4_mse-l1_pcc/log/339-formal.json) ·
[Cosine acceptance](05.5_mse-l1_cos/log/339-formal.json)
