# Step 06 — PerturbDiff transfer diagnostics

Step 05 observed a Tahoe Transfer deficit against the tenv-only CM arm.
Step 06 establishes a shared validation baseline and tests positive-dose
initialization and target base learning rate in separate paired arms.
[ADR 0048](../../docs/adr/0048-perturbdiff-transfer-diagnosis.md)
records the screening, confirmation and attribution rules.

## Validation baseline

Use the arc Python interpreter for verification and scoring. Use the pinned
PerturbDiff interpreter for prediction. Submit these phases to the dedicated
scheduler. The baseline reads sealed predecessor artifacts in place.

```bash
step=stage2/06_pdiff_transfer_optimize
cli="$step/execution/script/baseline.py"
pdiff=stage2/_envs/perturbdiff/bin/python
out="$step/output/baseline/fold-0-seed-17"

python "$cli" verify "$out" --seed 17
plan=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$cli" predict "$out" --seed 17 --plan-sha256 "$plan"
python "$cli" score "$out" --seed 17 --plan-sha256 "$plan"
```

The plan pins the selected EMA states and the common target-free validation
request. Both arms use the same validation well and group order, same-plate
controls, released Morgan fingerprints, and inference seed and batch layout.
The score phase joins frozen validation Δ after verifying both prediction
seals. `scores/` contains per-group directions and undefined reasons,
per-dose and pooled summaries, and the paired CM minus Transfer gaps. The
same commands accept seeds 42 and 2026 with their own output directories.
Use `--runtime /absolute/path/to/arc` to locate predecessor artifacts and the
pinned interpreter from a separate checkout.
