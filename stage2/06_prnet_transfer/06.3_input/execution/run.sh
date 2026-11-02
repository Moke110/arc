#!/usr/bin/env bash
# Stage 2 · 06.3 — build/verify the frozen two-source 18,164-gene input.
#
#   execution/run.sh --check                  # frozen sources + handoffs, no writes
#   execution/run.sh --out-root <step-dir>    # build + freeze 06.3_input/output/
#   execution/run.sh --verify --out-root <step-dir>
#
# Interpreter: the `arc` env (anndata 0.12.16) is the *writer*; the pinned
# trainer env (anndata 0.10.8) is the reader and cannot write the obs encoding.
# Override with PY=/path/to/python if the env moves.
set -euo pipefail

ROOT=$(pwd)
while [ "$ROOT" != "/" ] && [ ! -d "$ROOT/.git" ]; do ROOT=$(dirname "$ROOT"); done
if [ ! -d "$ROOT/.git" ]; then
  echo "error: not inside the ARC repository (no .git at or above $PWD)" >&2
  exit 1
fi
cd "$ROOT"

PY=${PY:-/home/chang/miniconda3/envs/arc/bin/python}
if [ ! -x "$PY" ]; then
  echo "error: python not found at $PY" >&2
  exit 1
fi

exec "$PY" stage2/06_prnet_transfer/06.3_input/execution/build_input.py "$@"
