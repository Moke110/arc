#!/usr/bin/env bash
# Stage 2 benchmark entry point.
#
#   bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --dry-run
#   bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model no_change --fold 0
#   bash stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/run.sh --model prnet --fold 0 --fold 1
#
# `--dry-run` runs the five CPU baselines on fold 0 and writes
# `stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/output/results/dry_run/leaderboard.csv`.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Repository root: walk up until a directory containing `.git` is found (ADR 0030 §4).
# Never count `../..` levels — a move must not silently change the root.
find_repo_root() {
  local dir="$1"
  while [ -n "$dir" ] && [ "$dir" != "/" ]; do
    if [ -e "$dir/.git" ]; then
      printf '%s\n' "$dir"
      return 0
    fi
    dir="$(dirname "$dir")"
  done
  echo "ERROR: no .git directory at or above $1" >&2
  return 1
}
ROOT="$(find_repo_root "$HERE")"
PY="${PYTHON:-/home/chang/miniconda3/envs/arc/bin/python}"

cd "$ROOT"
exec "$PY" "$HERE/runner.py" "$@"
