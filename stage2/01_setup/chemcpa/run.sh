#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OUT="$HERE/output"
ENV_NAME="${MODEL_ENV:-arc-stage2-bringup-chemcpa}"

mkdir -p "$OUT"
rm -f "$OUT/blocked.txt"
cd "$ROOT"
conda run --no-capture-output -n "$ENV_NAME" \
  python "$HERE/smoke.py" 2>&1 | tee "$OUT/run.log"

test -s "$OUT/pred.npy"
echo "chemCPA smoke OK: $OUT/pred.npy"
