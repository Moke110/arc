#!/usr/bin/env bash
# Sample total GPU memory during a benchmark queue (ticket-16 provenance).
# Writes `timestamp_s,total_used_mb` every 30 s; `summarize_run.py` attributes
# the peak to each unit by its status.json start/finish window.
set -euo pipefail
OUT="${1:-/home/chang/arc/stage2/03_bench_models/output/gpu_watch.csv}"
echo "timestamp_s,total_used_mb" > "$OUT"
while true; do
  ts=$(date +%s)
  mb=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1 | tr -d ' ')
  echo "${ts},${mb}" >> "$OUT"
  sleep 30
done
