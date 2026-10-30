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
PY="$ROOT/stage2/_envs/cmonge/bin/python"
mkdir -p "$OUT"
if [ ! -x "$PY" ]; then
  echo "missing $PY; see pert_model/CMonge/author_log.json" >&2
  exit 2
fi
"$PY" "$HERE/smoke.py" 2>&1 | tee "$OUT/run.log"
test -s "$OUT/meta.json"
echo "CMonge setup smoke OK: $OUT/meta.json"
