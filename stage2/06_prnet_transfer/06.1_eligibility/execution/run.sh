#!/usr/bin/env bash
# Stage 2 · 06 eligibility — reproduce the 460-structure eligibility audit.
#
#   bash execution/run.sh            # fetch/refresh evidence, build, run the gates
#   bash execution/run.sh --offline  # rebuild from the local evidence cache only
#   bash execution/run.sh --check    # gates on the committed record, no network
#
# Interpreter: the shared `arc` env (pandas + rdkit); no GPU is used.
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

exec "$PY" stage2/06_prnet_transfer/06.1_eligibility/execution/build_eligibility.py "$@"
