#!/usr/bin/env bash
# Sealed history is checked through its release record, separate from public routing.
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
cd "$root"
status=0
check_retired_paths() {
  local search_status=0
  rg -n "$@" || search_status=$?
  if [[ $search_status -gt 1 ]]; then
    echo 'ERROR: CPSP path search failed; no audit verdict is valid' >&2
    exit 2
  fi
  if [[ $search_status == 0 ]]; then status=1; fi
}
for path in cpsp/AGENTS.md cpsp/README.md cpsp/01_corpus/AGENTS.md cpsp/01_corpus/README.md cpsp/02_base_morgan/AGENTS.md cpsp/02_base_morgan/README.md pheno_data/AGENTS.md pheno_data/README.md docs/roadmap.md docs/data-inventory.md; do
  check_retired_paths 'stage3/|cpsp/01_data_corpus|\x6001_data_corpus/' "$path"
done
check_retired_paths 'stage3/|cpsp/01_data_corpus' cpsp/01_corpus/execution cpsp/02_base_morgan/execution --glob '*.py' --glob '!migrate_release.py'
if [[ -e stage3 || -L stage3 || -d cpsp/01_data_corpus || ! -f cpsp/01_corpus/execution/corpus.py || ! -f cpsp/02_base_morgan/execution/workflow.py ]]; then
  echo 'ERROR: CPSP must use its current paths' >&2
  status=1
fi
if [[ $status == 0 ]]; then echo 'CPSP source-first path audit: PASS'; fi
exit "$status"
