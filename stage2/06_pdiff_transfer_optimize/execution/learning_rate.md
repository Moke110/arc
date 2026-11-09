# Base-learning-rate intervention

Run the seed 17 screen on cold-drug fold 0 validation. Each pair uses the
original Tahoe EMA, bridge, recipient and native recipe. The intervention
changes target AdamW base learning rate from 0.0002 to 0.0001. Checkpoints
from the intervention are numbered under `output/checkpoints/`.

Use the [baseline commands](../README.md#validation-baseline) to seal the
matching seed's validation plan and both original predictions first. Use
the arc Python interpreter for verification, scoring and the decision. Use
the pinned PerturbDiff interpreter for fitting and prediction. The scheduler
admits each phase under the repository resource limits.

```bash
step=stage2/06_pdiff_transfer_optimize
cli="$step/execution/script/learning_rate.py"
pdiff=stage2/_envs/perturbdiff/bin/python
out="$step/output/learning-rate"
seed=17
baseline="$step/output/baseline/fold-0-seed-$seed/plan.json"
baseline_sha=$(sha256sum "$baseline" | cut -d ' ' -f 1)
python "$cli" verify "$out" --seed "$seed" \
  --baseline-plan "$baseline" --baseline-plan-sha256 "$baseline_sha"
plan=$(sha256sum "$out/fold-0-seed-$seed/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$cli" fit "$out" --seed "$seed" --plan-sha256 "$plan"
"$pdiff" "$cli" predict "$out" --seed "$seed" --plan-sha256 "$plan"
python "$cli" score "$out" --seed "$seed" --plan-sha256 "$plan"
python "$cli" decision "$out"
```

The plan predeclares the common optimizer steps based on the original
Transfer checkpoints. If seed 17 gains at least 0.010 overall validation
directional mean AP@100 at step 6500, repeat the commands for seed 42 at
step 6000 and seed 2026 at step 6500 using their own sealed baseline plans.
Reuse the original Transfer EMA at the declared step for each seed. Both
added seeds must gain strictly positive AP@100 to promote. The seed 17 CM
EMA matches its common step. The other seeds report CM under native selection
and mark equal-step CM recovery unavailable.

`native-selection.json` records the stop and the loss-selected EMA digest.
If native stopping occurs before the common step, the trainer seals the
native selection and continues the same optimizer, scheduler, random
stream and EMA to the common step under the original wall-time limit. The
diagnostic observations stay out of the native selection history.
`diagnostic.json` identifies the fixed-step EMA. `scores/` records both
comparison modes, directional and pooled AP@100 by dose, and coverage.
The terminal decision uses only the common-step pooled gain. If an admitted
seed ends before it can produce a score seal, declare its terminal reason
with `decision --unavailable SEED 'reason'`. A seed with a score seal uses
the scored common-step disposition.

The CPU-native fixture compares both native selection scenarios against the
untouched trainer and checks exact EMA and reload parity. The scheduler
admitted the fixture at two threads, with 10 GiB memory and 8 GiB swap limits.
It finished in 28.39 seconds and peaked at 4,182,272 KiB RSS.

```bash
"$pdiff" "$step/execution/tests/native_lr_fixture.py" --runtime "$PWD"
```
