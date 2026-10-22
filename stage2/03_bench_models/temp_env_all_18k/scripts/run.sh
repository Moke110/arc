#!/usr/bin/env bash
# Stage 2 benchmark entry point.
#
#   bash stage2/03_bench_models/temp_env_all_18k/scripts/run.sh --dry-run
#   bash stage2/03_bench_models/temp_env_all_18k/scripts/run.sh --model no_change --fold 0
#   bash stage2/03_bench_models/temp_env_all_18k/scripts/run.sh --model prnet --fold 0 --fold 1
#
# `--dry-run` runs the five CPU baselines on fold 0 and writes
# `stage2/03_bench_models/temp_env_all_18k/results/dry_run/leaderboard.csv`.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../../.." && pwd)"
PY="${PYTHON:-/home/chang/miniconda3/envs/arc/bin/python}"

cd "$ROOT"
exec "$PY" "$HERE/runner.py" "$@"
