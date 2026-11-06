#!/usr/bin/env bash
# User-authorized partial snapshot; source run tree is read-only. Run from repo root.
# A fresh output directory is required; never replace/forge the final 63-run seal.
set -euo pipefail
free -g
if [ "$(awk '/MemAvailable:/{print int($2/1048576)}' /proc/meminfo)" -lt 12 ]; then
  echo 'Insufficient available RAM' >&2
  exit 1
fi
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- \
  env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
  /usr/bin/time -v /home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/partial_benchmark.py \
  --input stage2/03_modes_tenv/output/input \
  --source stage2/03_modes_tenv/output/formal \
  --output stage2/03_modes_tenv/output/partial_finished/20261106_snapshot \
  --expected stage2/03_modes_tenv/execution/formal/expected_identities.json \
  --commands stage2/03_modes_tenv/execution/formal/commands.json --workers 2
# Lightweight checked rendering/publication; replaces the fixed visible figure paths.
env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
  /usr/bin/time -v /home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/partial_deliver.py \
  --snapshot stage2/03_modes_tenv/output/partial_finished/20261106_snapshot \
  --publish stage2/03_modes_tenv/output
