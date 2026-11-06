#!/usr/bin/env bash
# Live Stage 3 code/routing must use the source-first step; historical evidence is immutable.
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
cd "$root"
status=0
for path in stage3/AGENTS.md stage3/README.md stage3/01_corpus/AGENTS.md stage3/01_corpus/README.md pheno_data/AGENTS.md pheno_data/README.md docs/roadmap.md docs/data-inventory.md; do
  if rg -n 'stage3/01_data_corpus|`01_data_corpus/' "$path"; then status=1; fi
done
if rg -n 'stage3/01_data_corpus' stage3/01_corpus/execution --glob '*.py'; then status=1; fi
if [[ -d stage3/01_data_corpus || ! -f stage3/01_corpus/execution/corpus.py ]]; then
  echo 'ERROR: Stage 3 live step must be 01_corpus' >&2
  status=1
fi
if [[ $status == 0 ]]; then echo 'stage3 source-first path audit: PASS'; fi
exit "$status"
