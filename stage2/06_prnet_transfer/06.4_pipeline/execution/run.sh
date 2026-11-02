#!/usr/bin/env bash
# Stage 2 · 06.4 — one-command entry point for the PRnet formal-transfer pipeline.
#
#   execution/run.sh train --stage finetune --fold 0 --out <run-dir> --init random
#   execution/run.sh predict --stage finetune --fold 0 --role test \
#       --checkpoint <run-dir>/checkpoint/model.pt --out <pred-dir>
#   execution/run.sh evaluate --stage finetune --fold 0 --role test \
#       --pred-dir <pred-dir> --out <eval-dir>
#   execution/run.sh run_queue --queue <units.json> --out-root <out> --workers 2
#   execution/run.sh make_model_card --run-dir <run-dir> [--eval-dir <eval-dir>]
#   execution/run.sh smoke            # Phase-B real-data smoke (needs the frozen input)
#   execution/run.sh fixture          # synthetic fixture for CPU tests
#
# Interpreter: the pinned GPU env (torch 2.0.1+cu118).  Tests use the `arc` env
# through the Makefile.  Override with PY=/path/to/python.
set -euo pipefail

ROOT=$(pwd)
while [ "$ROOT" != "/" ] && [ ! -d "$ROOT/.git" ]; do ROOT=$(dirname "$ROOT"); done
if [ ! -d "$ROOT/.git" ]; then
  echo "error: not inside the ARC repository (no .git at or above $PWD)" >&2
  exit 1
fi

PY=${PY:-/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python}
if [ ! -x "$PY" ]; then
  echo "error: python not found at $PY" >&2
  exit 1
fi

SUB=stage2/06_prnet_transfer/06.4_pipeline/execution
cmd="${1:-}"
shift || true
case "$cmd" in
  train|predict|evaluate|run_queue|make_model_card)
    exec "$PY" "$ROOT/$SUB/$cmd.py" "$@"
    ;;
  smoke)
    exec "$PY" "$ROOT/$SUB/smoke_real.py" "$@"
    ;;
  fixture)
    exec "$PY" "$ROOT/$SUB/tools/make_fixture.py" "$@"
    ;;
  ""|-h|--help|help)
    sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
    ;;
  *)
    echo "error: unknown command '$cmd' (train|predict|evaluate|run_queue|make_model_card|smoke|fixture)" >&2
    exit 2
    ;;
esac
