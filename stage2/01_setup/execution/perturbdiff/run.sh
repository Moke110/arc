#!/usr/bin/env bash
set -euo pipefail
cursor="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
while [[ ! -e "$cursor/.git" ]]; do
  cursor="$(dirname -- "$cursor")"
done
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
"$cursor/stage2/_envs/perturbdiff/bin/python" "$cursor/stage2/01_setup/execution/perturbdiff/smoke.py" 2>&1 | tee "$cursor/stage2/01_setup/execution/perturbdiff/output/run.log"
