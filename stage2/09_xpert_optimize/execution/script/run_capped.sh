#!/usr/bin/env bash
# Run a Step 09 job serially with a bounded WSL scope and resource evidence.
set -euo pipefail

if (( $# < 2 )); then
    printf 'Usage: %s RECORD_DIR COMMAND [ARGS...]\n' "$0" >&2
    exit 2
fi
record=$1
shift
mkdir -p "$(dirname -- "$record")"
mkdir -- "$record"
exec 9>"${XDG_RUNTIME_DIR:-/run/user/$(id -u)}/arc-step09-heavy.lock"
if ! flock -n 9; then
    printf 'A Step 09 job already holds the heavy-job lock\n' | tee "$record/status.log" >&2
    exit 1
fi

available_kib() {
    awk '$1 == "MemAvailable:" {print $2}' /proc/meminfo
}

minimum_kib=$((12 * 1024 * 1024))
free -g | tee "$record/free-before.log"
if (( $(available_kib) < minimum_kib )); then
    printf 'Preflight failed. Available memory is below 12 GiB.\n' | tee "$record/status.log" >&2
    exit 1
fi

printf '%q ' "$@" > "$record/command"
printf '\n' >> "$record/command"
scope="arc-step09-$$.scope"
printf 'timestamp\tavailable_kib\tgpu_used_mib_total_mib_utilization_pct\tscope_peak_bytes\tworker_peak_rss_kib\n' \
    > "$record/resources.tsv"
(
    cgroup=""
    while :; do
        available=$(available_kib)
        gpu=$(timeout 2 nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu \
            --format=csv,noheader,nounits 2>> "$record/monitor-errors.log" || true)
        if [[ -z "$cgroup" ]]; then
            cgroup=$(systemctl --user show "$scope" --property=ControlGroup --value 2>/dev/null || true)
        fi
        scope_peak=""
        worker_peak=""
        if [[ -n "$cgroup" && -d "/sys/fs/cgroup$cgroup" ]]; then
            scope_peak=$(cat "/sys/fs/cgroup$cgroup/memory.peak" 2>/dev/null || true)
            worker_peak=$(
                while read -r pid; do
                    awk '/^VmHWM:/ {print $2}' "/proc/$pid/status" 2>/dev/null || true
                done < "/sys/fs/cgroup$cgroup/cgroup.procs" 2>/dev/null |
                    awk '$1 > maximum {maximum=$1} END {print maximum}'
            )
        fi
        printf '%s\t%s\t%s\t%s\t%s\n' "$(date -Is)" "$available" "$gpu" "$scope_peak" "$worker_peak" \
            >> "$record/resources.tsv"
        if (( available < minimum_kib )) && systemctl --user is-active --quiet "$scope"; then
            printf '%s Available RAM below 12 GiB. Interrupting the job.\n' "$(date -Is)" >> "$record/status.log"
            systemctl --user kill --signal=SIGINT "$scope" || true
            sleep 10
            if systemctl --user is-active --quiet "$scope"; then
                systemctl --user stop "$scope" || true
            fi
            break
        fi
        sleep 2
    done
) 9>&- &
monitor=$!
cleanup() {
    kill "$monitor" 2>/dev/null || true
    wait "$monitor" 2>/dev/null || true
}
trap cleanup EXIT

result=0
started=$EPOCHREALTIME
exec {stdout_fd}> >(tee "$record/stdout.log" 9>&-)
stdout_writer=$!
exec {stderr_fd}> >(tee "$record/stderr.log" >&2 9>&- {stdout_fd}>&-)
stderr_writer=$!
systemd-run --user --scope --unit="$scope" -p MemoryMax=10G -p MemorySwapMax=8G -- \
    /usr/bin/time -v -o "$record/time.log" env \
    OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 "$@" \
    9>&- 1>&"$stdout_fd" 2>&"$stderr_fd" {stdout_fd}>&- {stderr_fd}>&- || result=$?
finished=$EPOCHREALTIME
exec {stdout_fd}>&- {stderr_fd}>&-
if ! wait "$stdout_writer"; then
    (( result != 0 )) || result=1
fi
if ! wait "$stderr_writer"; then
    (( result != 0 )) || result=1
fi
printf '%s\t%s\n' "$started" "$finished" > "$record/wall-clock.tsv"
printf '%s exit=%d\n' "$(date -Is)" "$result" >> "$record/status.log"
free -g > "$record/free-after.log"
if (( result >= 128 )); then
    (dmesg 2>&1 || true) | tail -30 > "$record/kernel-tail.log"
fi
exit "$result"
