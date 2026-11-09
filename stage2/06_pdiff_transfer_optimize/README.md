# Step 06 — PerturbDiff transfer diagnostics

Step 05 observed a Tahoe Transfer deficit against the tenv-only CM arm.
Step 06 establishes a shared validation baseline and tests positive-dose
initialization and target base learning rate in separate paired arms.
[ADR 0048](../../docs/adr/0048-perturbdiff-transfer-diagnosis.md)
records the screening, confirmation and attribution rules.

## Results

Cold-drug fold 0 validation, seed 17. Each comparison covers 138 drug–dose
groups, with 46 at each of 0.2, 1 and 10 µM. Both AP@100 directions receive
equal weight. CM and original Transfer use the loss-selected step 6500 EMAs.

| Comparison | Intervention EMA step | Mean AP@100 | Gain vs original Transfer | Residual CM gap |
| --- | ---: | ---: | ---: | ---: |
| CM | 6500 | 0.25409 | | |
| Original Transfer | 6500 | 0.23843 | | 0.01566 |
| Positive-dose reset, loss-selected | 6000 | 0.24691 | +0.00848 | 0.00718 |
| Half base learning rate, matched step | 6500 | 0.23680 | −0.00163 | 0.01729 |
| Half base learning rate, loss-selected | 6000 | 0.23376 | −0.00467 | 0.02033 |

The dose reset improved the 0.2 and 1 µM strata by +0.01155 and +0.01470
respectively, with a −0.00080 change at 10 µM. Neither intervention reached
the +0.010 pooled promotion threshold. The added screening seeds and the
fold 1 source fits, target fits and test scoring were skipped. The dose
reset reduced the seed 17 validation gap by about 54%. The major-contribution
rule requires independent fold 1 replication and a matched gap comparison
there, so attribution remains unresolved.

Cross-dataset dose representation is worth further optimization. Replacing
the Tahoe-mapped positive-dose embeddings with fresh tenv embeddings
improved validation AP@100 in this screen, with gains concentrated at 0.2
and 1 µM. This single-seed result motivates further work on dose
representation and leaves its contribution to the Transfer deficit open.
The [fine-tuning research note](transfer_finetuning_research.md) compares
published adaptation methods with the original Transfer recipe.

The separate source audit covers nine shared tenv training parents and
50 Tahoe cell lines. Seven of those parents share identical standardized
parent SMILES. Historical Tahoe training overlapped one tenv validation
parent. After global exclusion, both Tahoe fitting roles overlap zero tenv
validation or test parents. Dose brackets and the assay, time, normalization
and coverage limitations are documented in the [source audit](execution/source_audit.md).

- [Paired effects figure](output/01-paired-validation-effects.png) and
  [Tahoe exposure figure](output/02-tahoe-exposure.png)
- [Paired per-run results](per_run.csv),
  [paired per-dose results](per_dose.csv) and
  [Tahoe fitting exposure](source_exposure.csv)
- [Decision summary](execution/publication_record/decision.json),
  [artifact hashes](execution/publication_record/evidence.json) and
  [resource ledger](execution/publication_record/resources.json)
- [Dose reset commands](execution/dose_reset.md),
  [learning-rate commands](execution/learning_rate.md),
  [confirmation source commands](execution/confirmation.md),
  [replication commands](execution/replication.md) and
  [publication replay](execution/publication.md)

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

## Paired Transfer reference

[ADR 0049](../../docs/adr/0049-perturbdiff-transfer-optimization-comparisons.md)
uses the Step 05 cold-cm-drug fold 0 original Transfer as a common reference
at seeds 17, 42 and 2026. Verify, predict and score each seed as separate jobs
through the dedicated scheduler. From the runtime checkout, use the arc
interpreter for verification and scoring and the pinned PerturbDiff
interpreter for prediction.

```bash
step=stage2/06_pdiff_transfer_optimize
cli="$step/execution/script/paired_reference.py"
pdiff="$PWD/stage2/_envs/perturbdiff/bin/python"
seed=17  # Repeat for 42 and 2026.
out="$step/output/paired-transfer/reference/seed-$seed"

python "$cli" verify "$out" --seed "$seed"
plan=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$cli" predict "$out" --seed "$seed" --plan-sha256 "$plan"
python "$cli" score "$out" --seed "$seed" --plan-sha256 "$plan"
```

The plan binds the original Tahoe and tenv EMAs and their author validation
loss histories. It records Tahoe's original exposure to tenv held-out drug
parents, along with the validation wells, group order, same-plate controls,
Morgan fingerprints, gene axis and sampling settings. The replay opens only
predictor inputs. The scoring job joins the sealed validation truth with the
sealed prediction and writes `scores/per_group.csv`, `scores/per_dose.csv`
and `scores/scores.json`. Dose summaries contain AUPRC (mean AP@100), PCC,
valid group counts and counts by undefined reason. The score record
distinguishes the selected EMA step from actual optimizer steps.

## Random positive-dose reset

The `random` arm restores the three native tenv positive-dose embeddings
after loading the same-seed Tahoe bridge. It then uses the original tenv
training recipe and selects its own EMA by the author validation diffusion
loss. Run all three seeds through the dedicated scheduler.

```bash
step=stage2/06_pdiff_transfer_optimize
cli="$step/execution/script/paired_target.py"
pdiff="$PWD/stage2/_envs/perturbdiff/bin/python"
seed=17  # Repeat for 42 and 2026.
reference="$step/output/paired-transfer/reference/seed-$seed"
out="$step/output/paired-transfer/random/seed-$seed"
reference_plan=$(sha256sum "$reference/plan.json" | cut -d ' ' -f 1)

python "$cli" verify "$out" --arm random --seed "$seed" --reference-sha256 "$reference_plan"
plan=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$cli" fit "$out" --arm random --seed "$seed" --plan-sha256 "$plan"
"$pdiff" "$cli" predict "$out" --arm random --seed "$seed" --plan-sha256 "$plan"
python "$cli" score "$out" --arm random --seed "$seed" --plan-sha256 "$plan"
```

The fit seal and selected checkpoint are under
`output/checkpoints/04-paired-random-fold-0-seed-$seed/fit/`.
The fit record includes the restored dose row digests, actual optimizer
steps and selected EMA step. The replay writes `prediction/`, the common
validation evaluator writes `scores/`, and `paired.json` records the
per-dose differences from original Transfer at the same seed.

## Continuous log dose

The `log` arm fits Tahoe on the original fold 0 source memberships and
recipe. The dose encoder applies a trainable affine map to
`log10(d / 0.05 µM)` and uses a separate trainable control vector at 0 µM.
Tenv starts from the same seed's loss-selected log Tahoe EMA and rebuilds
its own dose and Morgan lookups. The encoder evaluates 10 µM directly.

```bash
step=stage2/06_pdiff_transfer_optimize
source_cli="$step/execution/script/paired_log_source.py"
target_cli="$step/execution/script/paired_target.py"
pdiff="$PWD/stage2/_envs/perturbdiff/bin/python"
seed=17  # Repeat for 42 and 2026.
reference="$step/output/paired-transfer/reference/seed-$seed"
out="$step/output/paired-transfer/log/seed-$seed"
reference_plan=$(sha256sum "$reference/plan.json" | cut -d ' ' -f 1)

"$pdiff" "$source_cli" fit "$out" --seed "$seed" --reference "$reference" --plan-sha256 "$reference_plan"
python "$target_cli" verify "$out" --arm log --seed "$seed" --reference-sha256 "$reference_plan"
plan=$(sha256sum "$out/plan.json" | cut -d ' ' -f 1)
"$pdiff" "$target_cli" fit "$out" --arm log --seed "$seed" --plan-sha256 "$plan"
"$pdiff" "$target_cli" predict "$out" --arm log --seed "$seed" --plan-sha256 "$plan"
python "$target_cli" score "$out" --arm log --seed "$seed" --plan-sha256 "$plan"
```

The dedicated scheduler runs the jobs in order. Source and target fit
evidence lives under `output/checkpoints/08-paired-log-tahoe-fold-0-seed-$seed/fit/`
and `output/checkpoints/05-paired-log-fold-0-seed-$seed/fit/`, respectively.
