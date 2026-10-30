#!/usr/bin/env bash
set -euo pipefail

# Keep the vendored `source/` trees byte-identical to their downloads (ADR 0030):
# the smoke imports model code from them, and a stray .pyc would break the
# item's manifest. Exported so every child (`conda run`, the venv python) inherits.
export PYTHONDONTWRITEBYTECODE=1

HERE="$(cd "$(dirname "$0")" && pwd)"
# Repository root: walk up to the directory that contains .git (ADR 0030 §4).
find_repo_root() {
  local dir="$1"
  while [ -n "$dir" ] && [ "$dir" != "/" ]; do
    if [ -e "$dir/.git" ]; then printf '%s\n' "$dir"; return 0; fi
    dir="$(dirname "$dir")"
  done
  return 1
}
ROOT="$(find_repo_root "$HERE")" || {
  echo "error: cannot find the repository root above $HERE" >&2
  exit 2
}
OUT="$HERE/output"
PY="$HERE/.venv/bin/python"
if [ ! -x "$PY" ]; then
  echo "missing $PY; run the biolord setup command from README.md" >&2
  exit 2
fi

mkdir -p "$OUT"
cd "$ROOT"
PYTHONPATH="$ROOT/pert_model/biolord/source/src${PYTHONPATH:+:$PYTHONPATH}" \
  "$PY" "$HERE/smoke.py" 2>&1 | tee "$OUT/run.log"

test -s "$OUT/run.log"
echo "biolord smoke OK: $OUT/run.log"
