#!/usr/bin/env bash
# Stage 2 · 05 — one-command PRnet pretrain reproduction.
#
#   execution/run.sh                 # full 100-epoch run on the 04 pretrain view
#   execution/run.sh --check         # contract check only (paths, shape, split)
#   execution/run.sh --smoke --view stage2/05_prnet_pretrain/execution/_smoke/tahoe_cm_schema_smoke.h5ad
#
# Interpreter: the pinned GPU env, never the shared `arc` env (torch 2.12 removed
# ReduceLROnPlateau(verbose=...) and its scanpy/numba cache is broken).  Override with
# PY=/path/to/python if the env moves.
set -euo pipefail

ROOT=$(pwd)
while [ "$ROOT" != "/" ] && [ ! -d "$ROOT/.git" ]; do ROOT=$(dirname "$ROOT"); done
if [ ! -d "$ROOT/.git" ]; then
  echo "error: not inside the ARC repository (no .git at or above $PWD)" >&2
  exit 1
fi
cd "$ROOT"

PY=${PY:-/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python}
if [ ! -x "$PY" ]; then
  echo "error: python not found at $PY" >&2
  exit 1
fi

exec "$PY" stage2/05_prnet_pretrain/execution/train_prnet.py "$@"
