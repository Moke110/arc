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
printf 'timestamp\tavailable_kib\tgpu_used_mib_total_mib_utilization_pct\n' > "$record/resources.tsv"
(
    while :; do
        available=$(available_kib)
        gpu=$(timeout 2 nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu \
            --format=csv,noheader,nounits 2>> "$record/monitor-errors.log" || true)
        printf '%s\t%s\t%s\n' "$(date -Is)" "$available" "$gpu" >> "$record/resources.tsv"
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
systemd-run --user --scope --unit="$scope" -p MemoryMax=10G -p MemorySwapMax=8G -- \
    /usr/bin/time -v -o "$record/time.log" env \
    OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 "$@" \
    9>&- > >(tee "$record/stdout.log") 2> >(tee "$record/stderr.log" >&2) || result=$?
printf '%s exit=%d\n' "$(date -Is)" "$result" >> "$record/status.log"
free -g > "$record/free-after.log"
if (( result >= 128 )); then
    (dmesg 2>&1 || true) | tail -30 > "$record/kernel-tail.log"
fi
exit "$result"
