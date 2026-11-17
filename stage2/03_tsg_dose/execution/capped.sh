#!/usr/bin/env bash
# capped.sh — run one heavy dose-input step under the repository resource caps.
#
# The documented form is a systemd user scope:
#   systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- CMD
# Some agent sessions have no reachable user bus ("Failed to connect to bus:
# No data available"). In that case the same 10 GiB ceiling is applied with
# RLIMIT_AS, threads are pinned to 2 and /usr/bin/time -v still reports the
# measured peak. Either way the caller passes the timed command explicitly:
#   execution/capped.sh /usr/bin/time -v python3 ...
set -euo pipefail
available_kib=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
if (( available_kib < 12 * 1024 * 1024 )); then
  echo "available memory below12GiB; refusing heavy command" >&2
  exit 75
fi
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
if systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- true >/dev/null 2>&1; then
  echo "caps: systemd-run user scope MemoryMax=10G MemorySwapMax=8G threads=2" >&2
  exec systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- "$@"
fi
echo "caps: systemd user bus unavailable; applying RLIMIT_AS=10GiB threads=2 locally" >&2
ulimit -v 10485760
exec "$@"
