# Positive-dose initialization screen

The reset restores only the three same-seed native tenv positive-dose
embedding rows after loading the Tahoe bridge. The zero-dose row and the
remaining recipient state are preserved. The scoring comparator is the
original Transfer validation EMA selected by native diffusion loss.

Run from the repository root after the shared validation baseline is sealed.
Submit verification, fitting, prediction and scoring to the Step 06 scheduler.

```bash
step=stage2/06_pdiff_transfer_optimize
cli="$step/execution/script/dose_reset.py"
pdiff=stage2/_envs/perturbdiff/bin/python
seed=17
base=$(sha256sum "$step/output/baseline/fold-0-seed-$seed/plan.json" | cut -d ' ' -f 1)
python "$cli" verify --seed "$seed" --baseline-plan-sha256 "$base"
plan=$(sha256sum "$step/output/dose_reset/fold-0-seed-$seed/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$cli" fit --seed "$seed" --plan-sha256 "$plan"
"$pdiff" "$cli" predict --seed "$seed" --plan-sha256 "$plan"
python "$cli" score --seed "$seed" --plan-sha256 "$plan"
python "$cli" decision
```

The score seal retains the paired effects, original and residual CM gaps,
both AP@100 directions and their group support. The screening decision
records the pooled promotion disposition and any authorized added seeds.
