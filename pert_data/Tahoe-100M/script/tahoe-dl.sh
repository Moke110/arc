#!/usr/bin/env bash
# tahoe-dl — start/stop the resumable Tahoe-100M full download as a systemd --user unit.
#
#   ./tahoe-dl.sh start     # start (or resume) in the background; returns immediately
#   ./tahoe-dl.sh pause     # SIGTERM: stop new work, drain in-flight files, exit
#   ./tahoe-dl.sh resume    # alias for start
#   ./tahoe-dl.sh status    # unit state + live status.json (files/bytes/rate/ETA)
#   ./tahoe-dl.sh logs [N]  # tail the download log
#
# The unit is transient (no unit file on disk) and survives this shell / DSH session.
# After a WSL restart, just run `start` again — the ledger makes it a no-op for done files.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ROOT="$SCRIPT_DIR"
while [ ! -e "$ROOT/.git" ]; do
  [ "$ROOT" = "/" ] && { echo "tahoe-dl: cannot find the repository root above $SCRIPT_DIR" >&2; exit 2; }
  ROOT="$(dirname -- "$ROOT")"
done
UNIT=tahoe-dl
STATUS="$ROOT/pert_data/Tahoe-100M/local/_state/status.json"

is_active() { systemctl --user is-active --quiet "$UNIT" 2>/dev/null; }

case "${1:-status}" in
  start|resume)
    if is_active; then
      echo "tahoe-dl: already running (systemctl --user is-active $UNIT)"
      exit 0
    fi
    systemd-run --user --unit="$UNIT" --collect \
      --working-directory="$ROOT" \
      --property=TimeoutStopSec=300 \
      --setenv="TAHOE_WORKERS=${TAHOE_WORKERS:-2}" \
      --description="Tahoe-100M full fetch (resumable)" \
      "$SCRIPT_DIR/run_fetch.sh"
    echo "tahoe-dl: started"
    ;;
  stop|pause)
    if ! is_active; then
      echo "tahoe-dl: not running"
      exit 0
    fi
    systemctl --user stop "$UNIT"
    echo "tahoe-dl: stopped (in-flight files drained; rerun start to resume)"
    ;;
  status)
    if is_active; then echo "unit: active"; else echo "unit: inactive"; fi
    echo "log:  $ROOT/logs/tahoe_full.log"
    if [ -f "$STATUS" ]; then
      echo "--- _state/status.json ---"
      cat "$STATUS"
      echo
    fi
    ;;
  logs)
    tail -n "${2:-40}" "$ROOT/logs/tahoe_full.log"
    ;;
  *)
    echo "usage: $0 start|pause|resume|stop|status|logs [N]" >&2
    exit 2
    ;;
esac
