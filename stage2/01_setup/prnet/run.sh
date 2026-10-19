#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OUT="$HERE/output"
ENV_NAME="${MODEL_ENV:-arc-stage2-bringup-prnet}"

mkdir -p "$OUT"
rm -f "$OUT"/demo_split_*_array.csv

cd "$ROOT/pert_model/prnet"
conda run --no-capture-output -n "$ENV_NAME" \
  python "$HERE/run_prnet.py" \
    --split_key demo_split \
    --data_path ./dataset/demo.h5ad \
    --save_dir ./checkpoint/ \
    --results_dir "$OUT/" 2>&1 | tee "$OUT/run.log"

test -s "$OUT/demo_split_y_pre_array.csv"
echo "PRnet smoke OK: $OUT/demo_split_y_pre_array.csv"
