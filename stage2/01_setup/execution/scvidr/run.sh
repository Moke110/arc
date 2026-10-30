#!/usr/bin/env bash
set -euo pipefail

# Keep the vendored `source/` trees byte-identical to their downloads (ADR 0030):
# the smoke imports model code from them, and a stray .pyc would break the
# item's manifest. Exported so every child (`conda run`, the venv python) inherits.
export PYTHONDONTWRITEBYTECODE=1

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/output"
ENV=arc-stage2-bringup-scvidr
mkdir -p "$OUT"
conda run -n "$ENV" python "$HERE/smoke.py" 2>&1 | tee "$OUT/run.log"
test -s "$OUT/meta.json"
test -s "$OUT/run.log"
echo "scVIDR setup smoke OK: $OUT/meta.json"
