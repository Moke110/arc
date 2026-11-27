# Step 06.1 TranSiGen expression shape

Read the [CarPS-Transcriptome routing](../../AGENTS.md), the [study report](README.md)
and [ADR 0094](../../../docs/adr/0094-transigen-shape-round1.md). For encoder-width
changes or reuse of prior fits, read [ADR 0095](../../../docs/adr/0095-transigen-shape-encoder-width-extension.md).
For runner or scoring changes, use the [shared TranSiGen APIs](../../util/README.md#transigen-training-and-evaluation).

## Register and fit

Use the frozen tenv cold-drug folds on the ordered Denoise-DEG 5k axis,
measured same-plate DMSO and seed 42. Materialize the [Step 03.7 training release](../../03_data_train/03.7_tenv_denoise-5k/output/manifest.json),
the frozen tenv KPGT vectors under `cpst/util/output/representations` and the
[read-only native TranSiGen source](../../../pert_model/TranSiGen/manifest.sha256).
For scoring, materialize [Step 03.1 HVG-DEG 5k](../../03_data_train/03.1_tenv_hvg-5k/output/manifest.json)
for fold bandwidth calibration. Preserve the objective, initialization,
optimizer and checkpoint selection in the [Step 05.2 weight-1.75 protocol](../../05_tsg_loss/05.2_mse_l1/README.md)
and Step 04 latent numeric dose conditioning. Apply each shape to both
expression branches with independent parameters. Keep the treated decoder
shared between reconstruction and prediction.

Use [`run_shape.py`](execution/run_shape.py) to register each new coordinate
before fitting. Run registration, training, prediction and independent
replay in separate processes. Apply the [execution environment](../../../docs/agents/execution-environment.md)
and define `capped_step04` in Bash using the [Step 04 runtime](../../04_tsg_axis/AGENTS.md#prepare-the-capped-runtime).
Run the commands below from the materialized repository root in `arc`.
Run one CUDA fit at a time and sample WSL RAM and GPU usage every ten seconds.
Qualify changes with the [bounded smoke](#reproduce-the-bounded-smoke).

## Verify the accepted #347 source

Keep the sealed matrix and its sibling `report` under `output/formal/prior-347/`.
Pass `--prior-matrix` if the matrix and report are materialized elsewhere.
Preserve the original coordinate names in registrations and receipts.
[Sealed source and names](README.md#verify-the-accepted-347-source)

```bash
python cpst/06_tsg_module/06.1_round1/execution/report_shape.py verify-prior
```

For a CUDA replay, use the registered environment and replay a scratch copy of
the fit. The example uses the #347 `widen_enc/fold-0` fit. For a sealed #350
extension fit, use `2a5ac0a7a6c7dba527820564455f638d7e834a81` and its run under
`output/formal/matrix/runs/`, for example `widen_enc-2k/fold-0`.
[#347 receipt](log/347-formal.json) · [#350 receipt](log/350-formal.json)

```bash
(
  set -e
  stage="$PWD/cpst"
  study=06_tsg_module/06.1_round1
  revision=1f166d2f0dee53444131eee1dd83481c8dedce4f
  run="$stage/$study/output/formal/prior-347/matrix/runs/widen_enc/fold-0"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-shape-prior.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  git worktree add --detach "$work/code" "$revision"
  trap 'git worktree remove "$work/code"' EXIT
  cp -aL --reflink=auto -- "$run" "$work/run"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib" PYTHONDONTWRITEBYTECODE=1
  capped_step04 python "$work/code/stage2/06_tsg_shape/06.1_round1/execution/run_shape.py" replay \
    --output "$work/run" \
    --training "$stage/03_data_train/03.7_tenv_denoise-5k/output" \
    --representations "$stage/util/output/representations" \
    --source "$PWD/pert_model/TranSiGen/source"
)
```

## Reproduce the bounded smoke

Run the original shape and the two added widths on fold 0 for two epochs.
Require independent replay of each fit, verification of all 27 fold scores
and byte-identical renders of the two score figures.
[Bounded evidence](README.md#reproduce-the-bounded-smoke)

```bash
(
  set -e
  step="$PWD/cpst/06_tsg_module/06.1_round1"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-shape-smoke.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib" "$work/runs"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib" PYTHONDONTWRITEBYTECODE=1
  shapes=(original widen_enc-2k widen_enc-4k)
  runs=()
  for shape in "${shapes[@]}"; do
    run="$work/runs/$shape"
    capped_step04 python "$step/execution/run_shape.py" register \
      --output "$run" --shape "$shape" --fold 0 --smoke
    runs+=(--run "$run")
  done
  for shape in "${shapes[@]}"; do
    for phase in train predict replay; do
      capped_step04 python "$step/execution/run_shape.py" "$phase" \
        --output "$work/runs/$shape"
    done
  done
  for phase in score verify; do
    capped_step04 python "$step/execution/evaluate_shape.py" "$phase" \
      --smoke "${runs[@]}" --output "$work/report" --scratch "$work/scratch"
  done
  for target in figures figures-replay; do
    python "$step/execution/render_shape.py" --smoke --data "$work/report" \
      --output "$work/$target"
  done
  diff -r "$work/figures" "$work/figures-replay"
)
```

## Reproduce the formal matrix

Use `fc680914193735519a3b150328e9240dc57b1f27` to verify the combined scientific
receipt and historical figures. That publication pins the renderer and retains
the `stage2/` code paths recorded by the sealed registrations. Materialize the
[accepted Step 05.2 weight-1.75 runs](../../05_tsg_loss/05.2_mse_l1/AGENTS.md#reproduce-the-original-formal-sweep)
for control parity. [Formal evidence](README.md#reproduce-the-formal-matrix) ·
[Original release recipe](https://github.com/Moke110/carps/blob/fc680914193735519a3b150328e9240dc57b1f27/stage2/06_tsg_shape/06.1_round1/README.md#reproduce-the-formal-matrix)

Set `mode=build` to create fresh registrations and run the six extension fits
with the pinned matrix driver. Reserve at least 8 GiB for the build and use a
backed volume under `tmp/` if needed. Keep report scratch and candidate publication
on one filesystem. Start the coordinator sampler in the same WSL distro with
access to `nvidia-smi`. For a resumed matrix, keep that sampler running across
attempts in the same WSL boot. Require exact replays, Step 05.2 control parity,
and verification of all fold scores and eligibility counts before release.

```bash
(
  set -e
  mode=verify
  stage="$PWD/cpst"
  study=06_tsg_module/06.1_round1
  formal="$stage/$study/output/formal"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "${TSG_SHAPE_WORK_ROOT:-$PWD/tmp}/tsg-shape-formal.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib" PYTHONDONTWRITEBYTECODE=1
  revision=fc680914193735519a3b150328e9240dc57b1f27
  git worktree add --detach "$work/code" "$revision"
  pinned="$work/code/stage2/06_tsg_shape/06.1_round1"
  sampler=
  cleanup() {
    if [ -n "$sampler" ]; then
      kill -TERM "$sampler" 2>/dev/null || true
      wait "$sampler" || true
    fi
    git worktree remove "$work/code"
  }
  trap cleanup EXIT

  training="$stage/03_data_train/03.7_tenv_denoise-5k/output"
  report="$work/report"
  if [ "$mode" = build ]; then
    matrix="$work/matrix"
    audit="$work/matrix-audit"
    publication="$work/publication"
    samples="$work/coordinator-resources.jsonl"
    python "$pinned/execution/sample_resources.py" --output "$samples" &
    sampler=$!
    while [ ! -s "$samples" ]; do
      kill -0 "$sampler" 2>/dev/null || exit 1
      sleep 1
    done
    capped_step04 python "$pinned/execution/matrix_shape.py" \
      --output "$matrix" --logs "$audit" --training "$training" \
      --representations "$stage/util/output/representations" \
      --source "$PWD/pert_model/TranSiGen/source"
    kill -TERM "$sampler"
    wait "$sampler"
    sampler=
  else
    matrix="$formal/matrix"
    audit="$formal/matrix-audit"
    publication="$pinned"
    cp -aL -- "$formal/report" "$report"
  fi

  args=(
    --matrix "$matrix" --matrix-audit "$audit"
    --prior-matrix "$formal/prior-347/matrix" --report "$report"
    --publication "$publication" --activity "$work/report-audit"
    --scratch "$work/scratch" --training "$training"
    --bandwidth-training "$stage/03_data_train/03.1_tenv_hvg-5k/output"
    --previous-step "$stage/05_tsg_loss/05.2_mse_l1"
  )
  if [ "$mode" = build ]; then
    capped_step04 python "$pinned/execution/report_shape.py" build "${args[@]}" \
      --coordinator-samples "$samples" --code-revision "$revision"
  fi
  capped_step04 python "$pinned/execution/report_shape.py" verify "${args[@]}"
)
```

## Replay the published figures

Use the renderer environment recorded in the [current figure receipt](log/figure-render.json).
Render from the three published plot CSVs and compare all PNGs byte for byte.
[Plot CSVs and figures](README.md#replay-the-published-figures)

```bash
(
  set -e
  step="$PWD/cpst/06_tsg_module/06.1_round1"
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-architecture-figures.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib" PYTHONDONTWRITEBYTECODE=1
  python "$step/figure/script/render_comparison.py" \
    --data "$step/figure/data" --output "$work/figures"
  for figure in learning-curve.png metrics.png metrics-paired.png; do
    cmp "$step/figure/$figure" "$work/figures/$figure"
  done
)
```
