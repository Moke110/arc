#!/usr/bin/env bash
# Stage 2 · 06.6 — one-command entry point for the ticket-#6 random-init CM matrix.
#
#   execution/run.sh queue      # write output/units.json (45 units)
#   execution/run.sh plan       # resume-aware dry run of the queue
#   execution/run.sh scoped     # THE heavy command: all 45 candidates under the
#                               # systemd-run 12 GiB scope, with 1 Hz host/GPU peaks
#   execution/run.sh record     # derive output/record/* from output/runs/*
#   execution/run.sh check      # re-derive + gate the record
#   execution/run.sh resume     # prove the resume semantics (45/45 skip, interrupted re-run)
#   execution/run.sh all        # scoped -> record -> check (the full delivery sequence)
#   execution/run.sh test       # offline unit tests (arc env)
#   execution/run.sh tables     # derived CSV table cache for the pinned env
#
# Every command also exists as a Makefile target; run.sh is the documented
# single entry point.  Override the interpreter/data dir/output root with
# `PY`, `GPU_PY`, `DATA_DIR`, `OUT`, `WORKERS`, `SCOPE`.
set -uo pipefail

ROOT=$(pwd)
while [ "$ROOT" != "/" ] && [ ! -d "$ROOT/.git" ]; do ROOT=$(dirname "$ROOT"); done
if [ ! -d "$ROOT/.git" ]; then
  echo "error: not inside the ARC repository (no .git at or above $PWD)" >&2
  exit 1
fi

SUB=stage2/06_prnet_transfer/06.6_cm_random
EXEC="$ROOT/$SUB/execution"
PY=${PY:-/home/chang/miniconda3/envs/arc/bin/python}
GPU_PY=${GPU_PY:-/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python}
DATA_DIR=${DATA_DIR:-$ROOT/stage2/06_prnet_transfer/06.3_input/output}
OUT=${OUT:-$ROOT/$SUB/output}
TABLE_CACHE=${TABLE_CACHE:-$OUT/input_cache}
WORKERS=${WORKERS:-2}
OMP_THREADS=${OMP_THREADS:-2}
SCOPE=${SCOPE:-systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=12G}
export PYTHONDONTWRITEBYTECODE=1

run_py() { "$PY" "$EXEC/$1" "${@:2}"; }
run_gpu() { "$GPU_PY" "$EXEC/$1" "${@:2}"; }

cmd="${1:-}"
shift || true
case "$cmd" in
  queue)
    exec "$PY" "$EXEC/build_queue.py" --data-dir "$DATA_DIR" --out-dir "$OUT" "$@"
    ;;
  plan)
    exec "$GPU_PY" "$EXEC/run_matrix.py" --data-dir "$DATA_DIR" --out-dir "$OUT" \
      --table-cache "$TABLE_CACHE" --workers "$WORKERS" --dry-run "$@"
    ;;
  scoped)
    read -r -a _scope <<<"$SCOPE"
    exec "${_scope[@]}" -- "$GPU_PY" "$EXEC/run_matrix.py" --data-dir "$DATA_DIR" \
      --out-dir "$OUT" --table-cache "$TABLE_CACHE" --workers "$WORKERS" \
      --omp-threads "$OMP_THREADS" --scope "$SCOPE" "$@"
    ;;
  record)
    exec "$PY" "$EXEC/write_record.py" --data-dir "$DATA_DIR" --out-dir "$OUT" "$@"
    ;;
  check)
    exec "$PY" "$EXEC/write_record.py" --check --data-dir "$DATA_DIR" --out-dir "$OUT" "$@"
    ;;
  resume)
    read -r -a _scope <<<"$SCOPE"
    exec "${_scope[@]}" -- "$GPU_PY" "$EXEC/run_matrix.py" --verify-resume \
      --data-dir "$DATA_DIR" --out-dir "$OUT" --table-cache "$TABLE_CACHE" "$@"
    ;;
  all)
    rc=0
    "$0" scoped "$@" || rc=$?
    "$0" record || rc=$?
    "$0" check || rc=$?
    exit "$rc"
    ;;
  tables)
    exec "$PY" "$ROOT/stage2/06_prnet_transfer/06.4_pipeline/execution/tools/convert_frozen_tables.py" \
      --data-dir "$DATA_DIR" --cache-dir "$TABLE_CACHE" "$@"
    ;;
  test)
    cd "$EXEC" && exec "$PY" -m unittest discover -s tests -v
    ;;
  ""|-h|--help|help)
    sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'
    ;;
  *)
    echo "error: unknown command '$cmd'" >&2
    exit 2
    ;;
esac
