# Fold 1 confirmation and Tahoe source

Both fold 0 screening decisions must be terminal. Freeze recomputes the
promotion gate from the paired score evidence. It retains the reason for
each nonpromoted arm. If both arms are nonpromoted, the frozen plan records
`skipped`, and verification checks that no confirmation source artifacts
exist. Preparation and fitting are disabled in that case.

Run from the repository root. Submit data verification and fitting to the
Step 06 scheduler. `freeze` uses `output/dose_reset/decision.json` and
`output/learning-rate/decision.json`. Use `--learning-rate-root` if the
learning-rate screen ran in a different directory.

```bash
runtime="$PWD"
step="$runtime/stage2/06_pdiff_transfer_optimize"
output="$step/output"
cli="$step/execution/script/confirmation.py"
pdiff="$runtime/stage2/_envs/perturbdiff/bin/python"

python "$cli" freeze "$output" --runtime "$runtime"
plan=$(sha256sum "$output/confirmation/plan.json" | cut -d ' ' -f 1)
python "$cli" verify "$output" --runtime "$runtime" --plan-sha256 "$plan"
```

For a ready plan, the confirmation recipient is the pinned Step 04 tenv
fold 1 release. Tahoe retains its Step 04 outer fold 0 and its per-seed
formal recipe. The union of the tenv fold 1 validation and test parent
connectivity blocks is excluded from both Tahoe fitting roles. The
Tahoe test role is retained. Preparation preserves the original case
order, control indices, gene axis, fingerprints, and expression rows.
It reads the Tahoe expression in bounded row batches. Before fitting and
independent replay, verification checks the filtered expression against
its ordered predecessor rows. It verifies both fitting roles separately.

```bash
python "$cli" prepare "$output" --runtime "$runtime" --plan-sha256 "$plan"
seed=17
"$pdiff" "$cli" fit "$output" --runtime "$runtime" --plan-sha256 "$plan" \
    --seed "$seed" --admission "$preflight"
"$pdiff" "$cli" reload "$output" --runtime "$runtime" --plan-sha256 "$plan" --seed "$seed"
python "$cli" verify-source "$output" --runtime "$runtime" --plan-sha256 "$plan" --seed "$seed"
```

Fit uses the Step 04 scheduler's admitted preflight record for the exact
source command, with the 18 GiB formal training cap. Use one job at a
time and retain scheduler telemetry. Repeat fit and reload for seeds 42
and 2026. Each seed has one native loss-selected source EMA shared by its
confirmation arms. `fit.json` records source membership, the native
validation selection, input and recipe hashes, admission, wall time,
and resource measurements. `source.json` binds the fit and separate
native replay of its exact validation loss and prediction probe.
The bridge may consume an EMA only after both seals verify. Source
checkpoint slots use ordinals 07–09 under `output/checkpoints/`.
Successful fits and reloads are reused on retry.

The exact equal-step learning-rate policy, if promoted, carries forward
optimizer steps 6500, 6000 and 6500 for seeds 17, 42 and 2026,
respectively. The dose reset, if promoted, uses native loss-selected EMA
checkpoints. The plan pins the target-free inference request and records
the test scoring metric before confirmation work begins.

Run the bounded release-to-reload fixture through the pinned interpreter.

```bash
"$pdiff" "$step/execution/tests/native_confirmation_fixture.py" --runtime "$runtime"
```
