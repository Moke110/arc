#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/output"
PY="$HERE/../../../pert_model/cmonge/.venv/bin/python"
mkdir -p "$OUT"
if [ ! -x "$PY" ]; then
  echo "missing $PY; see pert_model/cmonge/VENDORED.txt" >&2
  exit 2
fi
"$PY" "$HERE/smoke.py" 2>&1 | tee "$OUT/run.log"
test -s "$OUT/meta.json"
echo "CMonge setup smoke OK: $OUT/meta.json"
