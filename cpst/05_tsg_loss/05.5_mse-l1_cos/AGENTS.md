# Step 05.5 — TranSiGen Delta cosine sweep

Read the [study report](README.md), [Step 05 routing](../AGENTS.md) and
[ADR 0092](../../../docs/adr/0092-transigen-bounded-delta-similarity-sweeps.md).
Use the [shared TranSiGen APIs](../../util/README.md#transigen-training-and-evaluation)
for runner and scoring changes.

## Inputs and protocol

Fit on the frozen [tenv Denoise-DEG 5k release](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json)
with [cold-drug folds 0–2](../../03_data_train/output/frozen-folds/tenv/manifest.json),
seed 42, [frozen tenv KPGT vectors](../../util/README.md) and the pinned
[TranSiGen source](../../../pert_model/TranSiGen/manifest.sha256). Score with
the [frozen HVG-DEG 5k bandwidths](log/bandwidths.csv).

Use measured same-plate DMSO as input and as the common reference for both
Delta vectors. Keep Step 05.3's unit reconstruction weight and Delta L1
weight 1.75. Sweep `alpha` over 0, 0.1, 0.2, 0.3, 0.4 and 0.5. Compute cosine
across the ordered genes on uncentered Delta vectors, using epsilon `1e-8`
and zero similarity for a zero-norm vector. Preserve the Step 05.3 objective
exactly at `alpha = 0`. Pair held-out score differences with the fresh
cosine zero-weight fit in the same fold and seed.
[Shared similarity objective](../../util/training/similarity_objective.py)

## Prepare the runtime

Activate `arc` at the materialized checkout root and define the
[Step 04 capped runtime](../../04_tsg_axis/AGENTS.md#prepare-the-capped-runtime).
Follow the [execution environment](../../../docs/agents/execution-environment.md).
Keep accepted inputs read-only and candidates under `tmp/`. Run one CUDA
fit at a time and record WSL RAM and GPU usage every ten seconds.

## Reproduce the smoke

Register both fold-0 fits before training and use `--metric cos` in each
runner and evaluator phase. For relocated inputs, the runner accepts
`--training`, `--representations` and `--source`. The evaluator accepts
`--training` and `--bandwidth-training`. Require exact replay of the two-epoch
fits, byte-identical score verification and two identical CSV-only renders.
[Runner](../similarity/run_similarity.py) ·
[Evaluator](../similarity/evaluate_similarity.py) ·
[Smoke renderer](execution/render.py) · [Smoke receipt](log/338-smoke.json)

```bash
(
  set -e
  loss="$PWD/cpst/05_tsg_loss"
  shared="$loss/similarity"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-cos-smoke.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib" "$work/csv-only"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib" PYTHONDONTWRITEBYTECODE=1
  runs=()
  for alpha in 0 0.5; do
    run="$work/alpha-$alpha"
    capped_step04 python "$shared/run_similarity.py" register \
      --metric cos --alpha "$alpha" --fold 0 --smoke --output "$run"
    runs+=(--run "$run")
  done
  for alpha in 0 0.5; do
    for phase in train predict replay; do
      capped_step04 python "$shared/run_similarity.py" "$phase" \
        --metric cos --output "$work/alpha-$alpha"
    done
  done
  for phase in score verify; do
    capped_step04 python "$shared/evaluate_similarity.py" "$phase" --metric cos \
      --smoke "${runs[@]}" --output "$work/report" --scratch "$work/scratch"
  done
  cp "$work/report/"{delta-metrics,fold-scores,epoch-objectives}.csv "$work/csv-only/"
  for attempt in first second; do
    capped_step04 python "$loss/05.5_mse-l1_cos/execution/render.py" --smoke \
      --data "$work/csv-only" --output "$work/$attempt"
  done
  for name in delta-metrics.png delta-metrics-improvement.png learning-curves.png; do
    cmp "$work/first/$name" "$work/second/$name"
  done
)
```

## Reproduce the formal sweep

Build and verify a fresh candidate with the
[joint PCC and cosine procedure](../README.md#reproduce-the-original-formal-delta-similarity-sweeps).
Require exact independent replay for all 18 cosine fits and verify the frozen
bandwidths. Match the three fresh zero-weight fits to
[sealed Step 05.3](../05.3_mse-l1_mse-x/log/336-formal.json) at reconstruction
weight 1 on objective traces, selected and stopping epochs, validation and test
predictions, Delta arrays, fold scores and eligibility counts. The
[formal receipt](log/339-formal.json) pins the joint matrix, preflight and
cosine report.

## Replay a sealed fit

Use the registered environment and the original `stage2/` paths at revision
`6e9031714000534714304b012959a20992feb214`. Copy the installed run into scratch
and pass the current locations of the same pinned inputs. This example replays
the cosine zero-weight fit on fold 0.
[Registration verification](../../util/training/verification.md)

```bash
(
  set -e
  root="$PWD"
  stage="$root/cpst"
  mkdir -p "$root/tmp"
  work=$(mktemp -d "$root/tmp/tsg-cos-replay.XXXXXX")
  git worktree add --detach "$work/code" 6e9031714000534714304b012959a20992feb214
  trap 'git worktree remove "$work/code"' EXIT
  mkdir "$work/scratch" "$work/run"
  sealed="$stage/05_tsg_loss/05.5_mse-l1_cos/output/formal/matrix/runs/cos/alpha-0/fold-0"
  cp --reflink=auto --dereference -- "$sealed"/* "$work/run/"
  export TMPDIR="$work/scratch" PYTHONDONTWRITEBYTECODE=1
  capped_step04 python "$work/code/stage2/05_tsg_loss/similarity/run_similarity.py" replay \
    --metric cos --output "$work/run" \
    --training "$stage/03_data_train/03.7_tenv_denoise-5k/output" \
    --representations "$stage/util/output/representations" \
    --source "$root/pert_model/TranSiGen/source"
)
```

## Verify the current figures

In `arc`, copy the three plotted CSVs to scratch and render twice. Require
each result to match the published PNG byte for byte. Use the
[joint verifier](../README.md#verify-the-current-figures) to check the render
receipts for all five Step 05 studies.
[CSV-only renderer](figure/script/render_comparison.py) ·
[Render receipt](log/figure-render.json)

```bash
(
  set -e
  figure="$PWD/cpst/05_tsg_loss/05.5_mse-l1_cos/figure"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-cos-figures.XXXXXX")
  mkdir "$work/csv-only"
  cp "$figure/data/"{learning-curve,metrics,metrics-paired}.csv "$work/csv-only/"
  export TMPDIR="$work" PYTHONDONTWRITEBYTECODE=1
  for attempt in first second; do
    MPLCONFIGDIR="$work/mpl-$attempt" python "$figure/script/render_comparison.py" \
      --data "$work/csv-only" --output "$work/$attempt"
    for name in learning-curve.png metrics.png metrics-paired.png; do
      cmp "$figure/$name" "$work/$attempt/$name"
    done
  done
)
```
