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
ENV_NAME="${MODEL_ENV:-arc-stage2-bringup-prnet}"

mkdir -p "$OUT"
rm -f "$OUT"/demo_split_*_array.csv

cd "$ROOT/pert_model/PRnet/source"
conda run --no-capture-output -n "$ENV_NAME" \
  python "$HERE/run_prnet.py" \
    --split_key demo_split \
    --data_path ./dataset/demo.h5ad \
    --save_dir ./checkpoint/ \
    --results_dir "$OUT/" 2>&1 | tee "$OUT/run.log"

test -s "$OUT/run.log"
echo "PRnet smoke OK: $OUT/run.log"
