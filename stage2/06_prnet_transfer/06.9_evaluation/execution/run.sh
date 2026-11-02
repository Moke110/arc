#!/usr/bin/env bash
# Canonical entry point for ticket #9 — evaluate the 30 selected checkpoints.
#
#   ./run.sh units           verify + write the 30-unit manifest (fails closed until #7)
#   ./run.sh units-random    harness dry run on the random branch only
#   ./run.sh plan            print the scoring plan (no test row is read)
#   ./run.sh score-valid     predict+evaluate the VALIDATION role only (harness smoke)
#   ./run.sh score-test      THE test run; needs CONFIRM_FROZEN_SELECTION=1
#   ./run.sh record          aggregate into output/record
#   ./run.sh cards           30 model cards
#   ./run.sh report          results section + README block
#   ./run.sh check           re-derive the record, byte-compare
#   ./run.sh all             units -> score-test -> record -> cards -> report -> check
#
# Frozen input, #6/#7 run trees and #8's record live in the main checkout; the
# heavy run is one unit at a time. Prefer `make score-test-scoped` (systemd-run
# 8 GiB cap). Override T09_DATA_DIR / T09_TABLE_CACHE / T09_OUT.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$here" && git rev-parse --show-toplevel)"
PY="${PYTHON:-/home/chang/miniconda3/envs/arc/bin/python}"
GPU_PY="${GPU_PYTHON:-/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python}"
STEP="$ROOT/stage2/06_prnet_transfer/06.9_evaluation"
DATA_DIR="${T09_DATA_DIR:-/home/chang/arc/stage2/06_prnet_transfer/06.3_input/output}"
TABLE_CACHE="${T09_TABLE_CACHE:-/home/chang/arc/stage2/06_prnet_transfer/06.4_pipeline/output/input_cache}"
BASELINES_DIR="${T09_BASELINES_DIR:-$ROOT/stage2/06_prnet_transfer/06.8_baselines/output}"
RANDOM_RECORD="${T09_RANDOM_RECORD:-$ROOT/stage2/06_prnet_transfer/06.6_cm_random/output/record/selected_checkpoints.csv}"
TAHOE_RECORD="${T09_TAHOE_RECORD:-$ROOT/stage2/06_prnet_transfer/06.7_cm_tahoe/output/record/selected_checkpoints.csv}"
OUT="${T09_OUT:-$STEP/output}"
UNITS="$OUT/record/units.json"

if [[ $# -lt 1 ]]; then
  sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
fi

common_units=(--data-dir "$DATA_DIR" --baselines-dir "$BASELINES_DIR" \
  --random-record "$RANDOM_RECORD" --tahoe-record "$TAHOE_RECORD" --out "$UNITS")
common_score=(--units "$UNITS" --out "$OUT" --data-dir "$DATA_DIR" \
  --table-cache "$TABLE_CACHE" --hash-cache "$OUT/hash_cache.json" --gpu-py "$GPU_PY")
common_record=(--units "$UNITS" --out "$OUT" --data-dir "$DATA_DIR" \
  --baselines-dir "$BASELINES_DIR")

case "$1" in
  units)
    exec "$PY" "$here/build_units.py" "${common_units[@]}"
    ;;
  units-random)
    exec "$PY" "$here/build_units.py" "${common_units[@]}" --allow-incomplete
    ;;
  plan)
    exec "$PY" "$here/score.py" "${common_score[@]}" --phase test --plan --confirm-frozen-selection
    ;;
  score-valid)
    exec "$PY" "$here/score.py" "${common_score[@]}" --phase valid --limit 2
    ;;
  score-test)
    if [[ "${CONFIRM_FROZEN_SELECTION:-}" != "1" ]]; then
      echo "refusing: export CONFIRM_FROZEN_SELECTION=1 only after the Lead confirms" \
           "#6/#7 selections are final and #5's sweep is done" >&2
      exit 4
    fi
    exec "$PY" "$here/score.py" "${common_score[@]}" --phase test --confirm-frozen-selection
    ;;
  record)
    exec "$PY" "$here/write_record.py" "${common_record[@]}" --phase test
    ;;
  record-valid)
    exec "$PY" "$here/write_record.py" "${common_record[@]}" --phase valid --allow-subset --no-figures
    ;;
  cards)
    exec "$PY" "$here/make_cards.py" --units "$UNITS" --out "$OUT" --phase test
    ;;
  report)
    exec "$PY" "$here/report.py" --record-dir "$OUT/record" --write-readme
    ;;
  check)
    exec "$PY" "$here/write_record.py" "${common_record[@]}" --phase test --check
    ;;
  all)
    "$0" units
    CONFIRM_FROZEN_SELECTION=1 "$0" score-test
    "$0" record
    "$0" cards
    "$0" report
    "$0" check
    ;;
  *)
    echo "error: unknown command '$1'" >&2
    echo "       (units|units-random|plan|score-valid|score-test|record|record-valid|cards|report|check|all)" >&2
    exit 2
    ;;
esac
