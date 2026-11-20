#!/usr/bin/env bash
# Run one Step 04 phase under a capped systemd scope, recording its evidence.
set -euo pipefail
if [[ $# -lt 3 ]]; then
  printf 'usage: %s EVIDENCE_ROOT PHASE COMMAND [ARG ...]\n' "$0" >&2
  exit 64
fi
base=$1
phase=$2
shift 2
if [[ ! "$phase" =~ ^[a-z0-9-]+$ ]]; then
  printf 'invalid phase: %s\n' "$phase" >&2
  exit 64
fi
nvidia=$(command -v nvidia-smi || printf '%s\n' /usr/lib/wsl/lib/nvidia-smi)
mkdir -p "$base"
exec 9>"$base/gpu.lock"
flock -n 9
available_kib=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
(( available_kib >= 12 * 1024 * 1024 ))
available_gpu=$("$nvidia" --query-gpu=memory.free --format=csv,noheader,nounits)
(( available_gpu >= 4096 ))
attempt="$base/phases/$phase-$(date +%s%N)"
mkdir -p "$attempt"
printf '%s\n' "$attempt"
printf 'started=%s\nram_kib=%s\ngpu_mib=%s\n' "$(date -Iseconds)" "$available_kib" "$available_gpu" > "$attempt/preflight.txt"
printf '%q ' "$@" > "$attempt/command.txt"
printf '\n' >> "$attempt/command.txt"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
unit="arc-stage2-transfer-${attempt##*/}"
systemd-run --user --scope --unit "$unit" -p MemoryMax=10G -p MemorySwapMax=8G -- \
  /usr/bin/time -v -o "$attempt/time.txt" "$@" > "$attempt/stdout.txt" 2> "$attempt/stderr.txt" 9>&- &
job_pid=$!
status=0
(
  pressure=0
  for _ in 1 2 3 4 5; do
    if systemctl --user show "$unit.scope" -p Id -p InvocationID -p MemoryMax -p MemorySwapMax \
        -p MemoryCurrent -p MemoryPeak > "$attempt/scope.txt" 2>/dev/null \
        && awk -F= '$1 == "InvocationID" && length($2) == 32 { found = 1 } END { exit !found }' "$attempt/scope.txt"; then
      break
    fi
    sleep 0.1
  done
  while kill -0 "$job_pid" 2>/dev/null; do
    {
      date -Iseconds
      awk '/^(MemTotal|MemAvailable|SwapFree):/' /proc/meminfo
      "$nvidia" --query-gpu=memory.used,memory.free,utilization.gpu --format=csv,noheader,nounits
    } >> "$attempt/resources-10s.log"
    ram=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
    if (( ram < 12 * 1024 * 1024 )); then
      pressure=$((pressure + 1))
    else
      pressure=0
    fi
    if (( pressure >= 3 )); then
      systemctl --user stop "$unit.scope"
      exit 0
    fi
    sleep 10
  done
) 9>&- &
monitor_pid=$!
wait "$job_pid" || status=$?
kill "$monitor_pid" 2>/dev/null || true
wait "$monitor_pid" 2>/dev/null || true
printf 'ended=%s\nstatus=%s\n' "$(date -Iseconds)" "$status" > "$attempt/terminal.txt"
exit "$status"
