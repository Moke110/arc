#!/usr/bin/env bash
# capped.sh — run one heavy step-owned command under the repository caps.
#
# The documented form is a systemd user scope:
#   systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- CMD
# Some agent sessions have no reachable user bus ("Failed to connect to bus:
# No data available").  In that case the same 10 GiB ceiling is applied with
# RLIMIT_AS and threads are pinned to 2.  Either way the caller passes the
# timed command explicitly:
#   execution/capped.sh /usr/bin/time -v /home/chang/miniconda3/envs/arc/bin/python ...
set -euo pipefail
# Recheck shared VM headroom at the actual heavy-command boundary.
available_gib="$(free -g | awk '/^Mem:/ {print $7}')"
if ! [[ "$available_gib" =~ ^[0-9]+$ ]] || [ "$available_gib" -lt 12 ]; then
  echo "Refusing heavy command: require at least 12 GiB available RAM (free -g reports ${available_gib:-unknown}). Wait for other jobs and retry." >&2
  exit 75
fi
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
if systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- true >/dev/null 2>&1; then
  echo "caps: systemd-run user scope MemoryMax=10G MemorySwapMax=8G threads=2" >&2
  exec systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- "$@"
fi
echo "caps: systemd user bus unavailable; applying RLIMIT_AS=10GiB threads=2 locally" >&2
ulimit -v 10485760
exec "$@"
