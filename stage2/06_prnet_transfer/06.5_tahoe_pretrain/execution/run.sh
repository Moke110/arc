#!/usr/bin/env bash
# Stage 2 · 06.5 — one-command entry point for the seeded Tahoe pretrains (ticket #5).
#
#   execution/run.sh formal                 # the three frozen seeds, locked budget, 18G cap
#   execution/run.sh formal --seeds 20261021 # one seed (writes a *partial* formal record)
#   execution/run.sh profile                # short real-data profile (smoke + one full epoch)
#   execution/run.sh check-formal           # validate output/formal/tahoe_pretrain_record.json
#   execution/run.sh check-profile          # validate output/profile/profile_record.json
#   execution/run.sh plan                   # print the formal plan and exit (no data, no GPU)
#   execution/run.sh test                   # offline unit tests (arc env)
#
# The formal/profile commands run inside `systemd-run --user --scope` with a
# MemoryMax/MemorySwapMax cap by default (`--cap 12G`, `--cap 18G`, `--no-cap`).
# The wrapper itself refuses to start when the frozen input, the derived memmap
# cache, the frozen-table cache, the locked recipe or the memory floor disagree
# with the 06.3 contract.
#
# Interpreter: the pinned GPU env (torch 2.0.1+cu118).  Override with PY=...
set -euo pipefail

ROOT=$(pwd)
while [ "$ROOT" != "/" ] && [ ! -e "$ROOT/.git" ]; do ROOT=$(dirname "$ROOT"); done
if [ ! -e "$ROOT/.git" ]; then
  echo "error: not inside the ARC repository (no .git at or above $PWD)" >&2
  exit 1
fi

GPU_PY=${PY:-/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python}
ARC_PY=${ARC_PY:-/home/chang/miniconda3/envs/arc/bin/python}
SUB=stage2/06_prnet_transfer/06.5_tahoe_pretrain/execution
OUT=stage2/06_prnet_transfer/06.5_tahoe_pretrain/output
DATA_DIR=${DATA_DIR:-$ROOT/stage2/06_prnet_transfer/06.3_input/output}
CACHE=${CACHE:-$ROOT/stage2/06_prnet_transfer/06.4_pipeline/output/cache/tahoe_log1p_cp10k_v1.npy}
TABLE_CACHE=${TABLE_CACHE:-$ROOT/stage2/06_prnet_transfer/06.4_pipeline/output/input_cache}

cmd="${1:-}"
shift || true

# pull the optional cap flags out of the passthrough arguments
ARGS=()
CAP_OVERRIDE=""
for arg in "$@"; do
  case "$arg" in
    --no-cap) CAP_OVERRIDE=none ;;
    --cap=*) CAP_OVERRIDE="${arg#--cap=}" ;;
    *) ARGS+=("$arg") ;;
  esac
done

cap_for_mode() {
  case "$1" in
    profile) printf '12G' ;;
    formal) printf '18G' ;;
    *) printf 'none' ;;
  esac
}

case "$cmd" in
  profile|formal)
    CAP=${CAP_OVERRIDE:-$(cap_for_mode "$cmd")}
    BASE=("$GPU_PY" "$ROOT/$SUB/run_pretrain.py" --mode "$cmd" --data-dir "$DATA_DIR" \
          --cache "$CACHE" --table-cache "$TABLE_CACHE")
    if [ "${#ARGS[@]}" -gt 0 ]; then BASE+=("${ARGS[@]}"); fi
    if [ "$CAP" = "none" ]; then
      echo "06.5 $cmd (uncapped: assuming an outer scope already caps this process)" >&2
      exec "${BASE[@]}"
    fi
    echo "06.5 $cmd inside systemd-run --user --scope -p MemoryMax=$CAP -p MemorySwapMax=$CAP" >&2
    exec systemd-run --user --scope -p "MemoryMax=$CAP" -p "MemorySwapMax=$CAP" -- "${BASE[@]}"
    ;;
  plan)
    exec "$GPU_PY" "$ROOT/$SUB/run_pretrain.py" --mode formal --dry-run \
      --data-dir "$DATA_DIR" --cache "$CACHE" --table-cache "$TABLE_CACHE" "$@"
    ;;
  check-formal)
    exec "$ARC_PY" "$ROOT/$SUB/check_record.py" \
      --record "$ROOT/$OUT/formal/tahoe_pretrain_record.json" --expect-mode formal "$@"
    ;;
  check-profile)
    exec "$ARC_PY" "$ROOT/$SUB/check_record.py" \
      --record "$ROOT/$OUT/profile/profile_record.json" --expect-mode profile "$@"
    ;;
  test)
    exec "$ARC_PY" -m unittest discover -s "$ROOT/$SUB/tests" -v
    ;;
  ""|-h|--help|help)
    sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
    ;;
  *)
    echo "error: unknown command '$cmd' (formal|profile|plan|check-formal|check-profile|test)" >&2
    exit 2
    ;;
esac
