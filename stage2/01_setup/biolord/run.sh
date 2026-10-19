#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OUT="$HERE/output"
PY="$HERE/.venv/bin/python"
if [ ! -x "$PY" ]; then
  echo "missing $PY; run the biolord setup command from README.md" >&2
  exit 2
fi

mkdir -p "$OUT"
cd "$ROOT"
PYTHONPATH="$ROOT/pert_model/biolord/src${PYTHONPATH:+:$PYTHONPATH}" \
  "$PY" "$HERE/smoke.py" 2>&1 | tee "$OUT/run.log"

test -s "$OUT/pred_mean.npy"
echo "biolord smoke OK: $OUT/pred_mean.npy"
