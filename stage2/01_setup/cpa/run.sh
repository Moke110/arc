#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OUT="$HERE/output"
ENV_NAME="${MODEL_ENV:-arc-stage2-bringup-cpa}"

mkdir -p "$OUT"
cd "$ROOT"
conda run --no-capture-output -n "$ENV_NAME" \
  python "$HERE/smoke.py" 2>&1 | tee "$OUT/run.log"

test -s "$OUT/cpa_pred.npy"
echo "CPA smoke OK: $OUT/cpa_pred.npy"
