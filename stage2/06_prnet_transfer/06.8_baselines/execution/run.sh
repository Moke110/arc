#!/usr/bin/env bash
# Canonical entry point for ticket #8 — five-fold simple baselines.
#
#   ./run.sh check     verify the frozen 06.3 input, write nothing
#   ./run.sh run       fit + score every fold, write the tracked record
#   ./run.sh fixture   build the synthetic frozen fixture (test aid only)
#   ./run.sh dry-run   full record on the synthetic fixture (pipeline check, not evidence)
#   ./run.sh mem-check row-replicated fixture under /usr/bin/time -v (memory gate)
#
# The frozen input and the compute output live in the main checkout (a worktree
# is removed after merge), so the data root is explicit.  Override with
# T08_DATA_DIR / T08_OUT_ROOT / T08_TABLE_CACHE; copy the tracked record back
# into this step with `make record`; the five record files are tracked through a
# narrow `.gitignore` exception (`06.8_baselines/output/`).
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$here" && git rev-parse --show-toplevel)"
PY="${PYTHON:-/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python}"
DATA_DIR="${T08_DATA_DIR:-$ROOT/stage2/06_prnet_transfer/06.3_input/output}"
OUT_ROOT="${T08_OUT_ROOT:-$ROOT/stage2/06_prnet_transfer}"
TABLE_CACHE="${T08_TABLE_CACHE:-$ROOT/stage2/06_prnet_transfer/06.4_pipeline/output/input_cache}"
FIXTURE="${T08_FIXTURE:-/tmp/t08-baselines-fixture}"

if [[ $# -lt 1 ]]; then
  sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
fi

case "$1" in
  check)
    exec "$PY" "$here/baselines.py" --check --data-dir "$DATA_DIR" \
      --table-cache "$TABLE_CACHE"
    ;;
  run)
    exec "$PY" "$here/baselines.py" --all --data-dir "$DATA_DIR" \
      --out-root "$OUT_ROOT" --table-cache "$TABLE_CACHE"
    ;;
  fixture)
    exec "${PYTHON:-/home/chang/miniconda3/envs/arc/bin/python}" \
      "$ROOT/stage2/06_prnet_transfer/06.4_pipeline/execution/tools/make_fixture.py" \
      --out "$FIXTURE"
    ;;
  dry-run)
    exec "${PYTHON:-/home/chang/miniconda3/envs/arc/bin/python}" \
      "$here/baselines.py" --all --data-dir "$FIXTURE" \
      --out-root /tmp/t08-baselines-dry-run --force
    ;;
  mem-check)
    exec make -C "$here" mem-check
    ;;
  *)
    echo "error: unknown command '$1' (check|run|fixture|dry-run|mem-check)" >&2
    exit 2
    ;;
esac
