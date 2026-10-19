#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/output"
ENV=arc-stage2-bringup-decant
mkdir -p "$OUT"
conda run -n "$ENV" python "$HERE/smoke.py" 2>&1 | tee "$OUT/run.log"
test -s "$OUT/meta.json"
test -s "$OUT/pred_gene_dose1.npy"
echo "DECANT setup smoke OK: $OUT/pred_gene_dose1.npy"
