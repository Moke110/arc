#!/usr/bin/env bash
# Ticket-16 five-fold queue: seven real models + five CPU baselines.
#
# - One GPU unit at a time: the whole queue holds the shared flock.
# - 8 GB cgroup RAM cap + 2 GB swap cap around the runner and all adapters
#   (ticket 10 §6.1); the kernel OOM-kills only this scope on a runaway unit.
# - Per-fit/predict timeout is the runner default (6 h); the global 16 h/model
#   cap is not enforced here because every measured five-fold estimate is well
#   below it (the largest is CPA at ~5-10 h).
# - Completed units are skipped on a restart; `--force` reruns them.
#
# Usage:
#   bash stage2/03_bench_models/run_queue.sh
#   bash stage2/03_bench_models/run_queue.sh --force
#   bash stage2/03_bench_models/run_queue.sh --model prnet --model cpa
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
PY="${PYTHON:-/home/chang/miniconda3/envs/arc/bin/python}"
LOCK=/tmp/opencode/stage2_bench_gpu.lock

# Fast/author-faithful models first, the two long runs (CPA, Squidiff) last.
# A caller-supplied `--model` list replaces the default queue.
args=()
has_model=0
for arg in "$@"; do
  if [ "$arg" = "--model" ]; then
    has_model=1
    break
  fi
done
if [ "$has_model" -eq 0 ]; then
  # Baselines first (cheap, always available), then fast models, then the two
  # long runs.
  MODELS=(no_change delta_mean delta_mean_dose ridge_fp nearest_drug \
          prnet chemcpa biolord perturbnet cmonge cpa squidiff)
  for model in "${MODELS[@]}"; do
    args+=(--model "$model")
  done
fi

mkdir -p "$(dirname "$LOCK")"
cd "$ROOT"
exec flock "$LOCK" systemd-run --user --scope \
  -p MemoryMax=8G -p MemorySwapMax=2G -- \
  env PYTHONDONTWRITEBYTECODE=1 \
      CUBLAS_WORKSPACE_CONFIG=:4096:8 \
      PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
      OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 \
      OPENBLAS_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4 \
  "$PY" "$HERE/runner.py" "${args[@]}" "$@"
