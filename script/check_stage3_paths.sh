#!/usr/bin/env bash
# Live Stage 3 code/routing must use the source-first step; historical evidence is immutable.
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
cd "$root"
status=0
check_retired_paths() {
  local search_status=0
  rg -n "$@" || search_status=$?
  if [[ $search_status -gt 1 ]]; then
    echo 'ERROR: Stage 3 path search failed; no audit verdict is valid' >&2
    exit 2
  fi
  if [[ $search_status == 0 ]]; then status=1; fi
}
for path in stage3/AGENTS.md stage3/README.md stage3/01_corpus/AGENTS.md stage3/01_corpus/README.md pheno_data/AGENTS.md pheno_data/README.md docs/roadmap.md docs/data-inventory.md; do
  check_retired_paths 'stage3/01_data_corpus|`01_data_corpus/' "$path"
done
check_retired_paths 'stage3/01_data_corpus' stage3/01_corpus/execution --glob '*.py'
if [[ -d stage3/01_data_corpus || ! -f stage3/01_corpus/execution/corpus.py ]]; then
  echo 'ERROR: Stage 3 live step must be 01_corpus' >&2
  status=1
fi
if [[ $status == 0 ]]; then echo 'stage3 source-first path audit: PASS'; fi
exit "$status"
