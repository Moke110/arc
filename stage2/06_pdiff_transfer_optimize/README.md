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
