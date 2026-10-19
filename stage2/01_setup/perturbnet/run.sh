#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/output"
ENV=arc-stage2-bringup-perturbnet
mkdir -p "$OUT"
conda run -n "$ENV" python "$HERE/smoke.py" 2>&1 | tee "$OUT/run.log"
test -s "$OUT/meta.json"
test -s "$OUT/pred_latent_lincs.npy"
echo "PerturbNet setup smoke OK: $OUT/meta.json"
