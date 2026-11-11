#!/usr/bin/env bash
# Execute the frozen study serially, then stage independently replayed figures.
set -euo pipefail
repo=$(git -C "$(dirname -- "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)
cd "$repo"
step=stage2/09_xpert_optimize
python=/home/chang/miniconda3/envs/arc/bin/python
plan=$step/execution/rerun_record/plan.json
policy=$step/execution/rerun_record/policy.json
work=$step/execution/tmp/rerun/formal
staging=$step/output/round1-staging
"$python" "$step/execution/script/rerun_queue.py" --plan "$plan" --policy "$policy" --work "$work"
"$python" "$step/execution/script/rerun_queue.py" --plan "$plan" --policy "$policy" --work "$work" --evaluate
"$step/execution/script/run_capped.sh" "$work/publication-resources" "$python" \
    "$step/execution/script/rerun_publication.py" score-publish --plan "$plan" \
    --work-root "$work" --predictions-root "$work/predictions" --loss-root "$work/test_loss" \
    --output "$staging"
"$python" "$step/execution/script/rerun_publication.py" replay --output "$staging"
printf '24 fits, test prediction replays and staged figure replay complete; awaiting visual/final acceptance.\n'
