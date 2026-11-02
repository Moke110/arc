#!/usr/bin/env bash
# Canonical entry point for the per-observation QC audit (ticket #2).
#
#   ./run.sh --lock     derive + freeze the thresholds  -> output/qc_rules.json
#   ./run.sh --build    apply the frozen rules          -> per-observation tables
#   ./run.sh --verify   re-hash every artifact and re-check the headline counts
#
# The frozen inputs (pert_data/, stage2/02|04 output/) live in the main checkout, not
# in a worktree, so the data root is explicit. Override with QC_DATA_ROOT /
# QC_OUT_ROOT; the tracked record is copied back into this step with `make record`.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-/home/chang/miniconda3/envs/arc/bin/python}"
DATA_ROOT="${QC_DATA_ROOT:-/home/chang/arc}"
OUT_ROOT="${QC_OUT_ROOT:-$DATA_ROOT}"

if [[ $# -lt 1 ]]; then
  echo "usage: $0 {--lock|--build|--verify}" >&2
  exit 2
fi

exec "$PY" "$here/build_qc.py" "$1" --data-root "$DATA_ROOT" --out-root "$OUT_ROOT"
