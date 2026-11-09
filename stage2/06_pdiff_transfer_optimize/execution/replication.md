# Fold 1 paired test replication

The replication plan requires the frozen confirmation plan. With no promoted
intervention, freeze and finalize record a skipped disposition and reject any
fold 1 fit or score artifacts. For promotion, finish the shared Tahoe source
preparation, native fits, and independent reloads in
[confirmation.md](confirmation.md) first.

```bash
runtime="$PWD"
step="$runtime/stage2/06_pdiff_transfer_optimize"
output="$step/output"
cli="$step/execution/script/replication.py"
pdiff="$runtime/stage2/_envs/perturbdiff/bin/python"

confirmation=$(sha256sum "$output/confirmation/plan.json" | cut -d ' ' -f 1)
python "$cli" freeze "$output" --runtime "$runtime" --confirmation-sha256 "$confirmation"
plan=$(sha256sum "$output/replication/plan.json" | cut -d ' ' -f 1)
python "$cli" status "$output" --runtime "$runtime" --plan-sha256 "$plan"
```

Submit real-data verification, fits, prediction, and scoring to the Step 06
scheduler. A ready plan pins all three Tahoe EMAs, the matching Step 04 CM
selections, the frozen train and validation recipe, the target-free ordered
test request, evaluator file hashes, and the per-seed common step. It records
prior Step 03 and Step 04 scoring and measures overlap with fold 0 roles at
parent, drug, well and group level. The Step 05 fold 0 score is historical
context.

Fit the original Transfer arm once per seed. Fit each promoted intervention
for the same seed, using a separate admitted command for each arm. Replace
`$preflight` with the admission record for the exact command. New checkpoints
use slots 10–12 under `output/checkpoints/`. Each pair verifies the same
source and mapped bridge. Dose reset restores the three native positive-dose
rows and keeps the original recipe. Learning rate changes only the base
rate. When learning rate is promoted, the original and halved-rate fits
also capture the frozen equal-step EMA from the same live trajectory.
Native loss selection remains separately available.

```bash
seed=17
arm=Transfer  # or a promoted dose_reset or learning_rate
"$pdiff" "$cli" fit "$output" --runtime "$runtime" --plan-sha256 "$plan" \
    --seed "$seed" --arm "$arm" --admission "$preflight"
# Run prediction and scoring after fitting the authorized arms for this seed.
"$pdiff" "$cli" predict "$output" --runtime "$runtime" --plan-sha256 "$plan" --seed "$seed"
python "$cli" score "$output" --runtime "$runtime" --plan-sha256 "$plan" --seed "$seed"
```

Repeat for seeds 42 and 2026. Prediction receives only cases, same-plate
controls, gene order and the released fingerprints. Scoring verifies the
prediction seals before opening frozen test truth and directional labels.
Per-seed scores include native and equal-step contrasts where applicable,
AP@100 by direction and dose, and the underlying group rows. Equal-step
attribution requires a CM snapshot at the same optimizer step. The
loss-selected CM appears in the native comparison when the equal-step CM is
unavailable.

```bash
python "$cli" finalize "$output" --runtime "$runtime" --plan-sha256 "$plan"
```

Finalize also handles skipped plans. For a seed without a score seal, record
a terminal disposition with `--unavailable SEED STATUS 'reason'` using
`failed`, `incomplete` or `skipped`. The final report at
`output/replication/result.json` retains seed evidence, primary paired
effects, equal-weight means, matching CM coverage, and any residual gaps.
A promotion replicates when all three primary seed pairs are available and
their mean intervention gain is positive. Recovery uses the comparable seed
set and requires a positive mean original CM gap.
