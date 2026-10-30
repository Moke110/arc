#!/usr/bin/env bash
# Wrapper for the systemd --user unit `tahoe-dl`.
# Runs the resumable full-scope Tahoe-100M fetch; stdout/stderr go to logs/tahoe_full.log.
# Safe to kill/restart at any time: fetch_tahoe.py has its own ledger + .incomplete resume.
set -euo pipefail
cd /home/chang/arc
export HF_ENDPOINT="${HF_ENDPOINT:-https://hf-mirror.com}"
mkdir -p logs
# Worker count default 2 (operator choice, 2026-10-29). Measured on hf-mirror:
# 1 stream ~5.4 MB/s, 2 streams ~7.0 MB/s in-service, 4 streams ~2-4 MB/s
# (per-connection throttling). Override with TAHOE_WORKERS=1|4.
WORKERS="${TAHOE_WORKERS:-2}"
exec /home/chang/miniconda3/bin/python3 \
  stage2/04_tahoe_pretrain/script/fetch_tahoe.py \
  --scope full --workers "$WORKERS" \
  "$@" >> logs/tahoe_full.log 2>&1
