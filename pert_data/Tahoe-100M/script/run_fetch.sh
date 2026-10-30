#!/usr/bin/env bash
# Wrapper for the systemd --user unit `tahoe-dl`.
# Runs the resumable full-scope Tahoe-100M fetch; stdout/stderr go to logs/tahoe_full.log.
# Safe to kill/restart at any time: fetch_tahoe.py has its own ledger + .incomplete resume.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ROOT="$SCRIPT_DIR"
while [ ! -e "$ROOT/.git" ]; do
  [ "$ROOT" = "/" ] && { echo "run_fetch: cannot find the repository root above $SCRIPT_DIR" >&2; exit 2; }
  ROOT="$(dirname -- "$ROOT")"
done
cd "$ROOT"
export HF_ENDPOINT="${HF_ENDPOINT:-https://hf-mirror.com}"
mkdir -p logs
# Worker count default 2 (operator choice, 2026-10-29). Measured on hf-mirror:
# 1 stream ~5.4 MB/s, 2 streams ~7.0 MB/s in-service, 4 streams ~2-4 MB/s
# (per-connection throttling). Override with TAHOE_WORKERS=1|4.
WORKERS="${TAHOE_WORKERS:-2}"
# Prefer this host's interpreter; degrade to PATH, then fail with a clear message.
PY="${PYTHON:-/home/chang/miniconda3/bin/python3}"
if [ ! -x "$PY" ]; then
  PY="$(command -v python3 || true)"
fi
if [ -z "$PY" ]; then
  echo "$(basename "$0"): no usable python3 — set PYTHON=/path/to/python3 (see environment.yml)" >&2
  exit 1
fi
exec "$PY" \
  pert_data/Tahoe-100M/script/fetch_tahoe.py \
  --scope full --workers "$WORKERS" \
  "$@" >> logs/tahoe_full.log 2>&1
