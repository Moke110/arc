#!/usr/bin/env bash
# Stage 2 v1 whole-transcriptome CM benchmark — end-to-end.
#
#   bash run.sh [--adapters no-change mean-effect ridge-fp nearest-drug]
#
# Requires the shared `arc` conda env on PATH (conda activate arc). Model
# adapters additionally require their own pinned envs (ADR 0021).
set -euo pipefail
cd "$(dirname "$0")"

echo "== audit =="
python audit_entrants.py

echo "== 01_bench_data =="
python 01_bench_data/build_bench_data.py

echo "== 02_baselines (in-sample sanity) =="
python 02_baselines/gate.py

echo "== 03_stub (runner harness) =="
python 03_stub/gate.py

echo "== 04_evaluate =="
python 04_evaluate/evaluate.py "$@"
