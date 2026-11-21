# Development repeat reliability and model residuals (#262)

The frozen Step 03 validation partitions provide the treatment observations for this
audit. Use each run's own train split for the Mean_cell reference. The fitted model
is the fixed 32-dimensional latent embedding checkpoint restored at minimum treated
validation MSE (training seed 42, split seed 226). Each prediction and validation
identity file is linked to the sealed Step 03 fit manifest. G_5000 uses ADR 0077
axis-local normalization. The same row can be in more than one validation scenario
or fold. Tables keep those partitions separate. Validation scores describe models
on the same partition used for checkpoint selection.

## Repeats

Group by chemical parent, physical dose and cell background within each source,
scenario and fold, then pair treated observations within that group. Read the
original vehicle row identities from the freeze. Check vehicle role, source item,
plate, cell line, vehicle metadata when present, and unique original source row.
The original-reference estimands use the exact frozen matched control for each
treated observation. Different plate groups with disjoint original vehicle sets
provide cross-batch evidence; plate effects and control strategy are confounded
in that comparison.

For same-batch repeat pairs, split the original vehicles into two disjoint
equal-size pools by alternating source row order. Drop the last vehicle when the
count is odd. Compare using the first shared pool and the second shared pool,
averaged, with the mean over both disjoint assignments on the identical treated
pair. This keeps pool sample count fixed. The construction measures sensitivity
to reusing original vehicle observations within a batch. Intra-batch common
variation can remain even with disjoint vehicle samples.

| Validation scenario | Source | Natural shared pairs | Natural disjoint pairs | Matched sensitivity, shared minus disjoint Δ PCC |
|---|---|---:|---:|---:|
| Cold-drug | tenv | 873 | 0 | +0.1275 to +0.1559 |
| Cold-drug | Tahoe | 150 | 1,546 | +0.3337 to +0.4712 |
| Cold-drug-dose | tenv | 582 | 0 | +0.1275 to +0.1559 |
| Cold-drug-dose | Tahoe | 100 | 950 | +0.3337 to +0.3766 |

Counts sum fold-dose cells within each scenario. Score ranges use nonempty cells.
In tenv, each sensitivity cell represents 31 parents. The natural independent-batch
sensitivity cannot be estimated because each observed tenv condition's treatment
repeats share its original vehicle group. In Tahoe the matched comparison covers
one parent in each nonempty cell, with 50 distinct conditions per such cell.
Natural disjoint controls cover one to seven parents per cell. For identifying
missing support, use the dose-specific counts in [repeat coverage](repeat_coverage.csv)
and defined denominators in the [repeat estimates](repeat_summary.csv). Singleton
conditions are the principal Tahoe coverage loss.

Average eligible pairs within condition, conditions within parent, then parent
means equally within each fold-dose cell. Percentile 95% intervals use 2,000 draws
over parents (seed 262). Intervals are omitted for a single parent. Intervals
hold the vehicle split fixed and describe sampling among observed parents under
the frozen partitions. Repeat PCC is a reproducibility diagnostic; an upper bound
on prediction performance would require a measurement-error model and stronger
independence assumptions.

## Residuals

Use the original matched control for conventional Δ PCC. For the residual score,
subtract the train-only Mean_cell prediction from predicted and observed treated
expression before computing gene-wise PCC. This equals subtracting the row-specific
Mean_cell Δ prediction from both Δ vectors, so the shared control cancels.
Mean_cell pools treated training profiles and deduplicated original training
vehicles within cell line. Residual PCC may still include dose or other systematic
effects shared among parents.

| Validation scenario | Source | Original Δ PCC | Residual PCC |
|---|---|---:|---:|
| Cold-drug | tenv | 0.2253–0.3762 | 0.1647–0.2672 |
| Cold-drug | Tahoe | 0.4991–0.5970 | 0.3987–0.5063 |
| Cold-drug-dose | tenv | 0.2283–0.3517 | 0.1402–0.2474 |
| Cold-drug-dose | Tahoe | 0.5000–0.5903 | 0.3969–0.4621 |

Means are profile-equal within parent, then parent-equal within each fold-dose
cell. All cells contain 31 tenv or 24 Tahoe parents. See [residual estimates](residual_summary.csv)
for profile counts, paired differences and parent-bootstrap intervals. Centering
changes both vectors and changes the reference estimand, so the score differences
combine control removal and Mean_cell adjustment.

## Replay

Use the accepted Step 03 pair and formal matrix as described in the step README.
From the repository root in the arc environment, with a fresh output directory:

```bash
STEP=stage2/03_tsg_dose
"$STEP/execution/capped.sh" /usr/bin/time -v python \
  "$STEP/execution/axis_local/control_reliability.py" \
  --inputs "$STEP/execution/checkpoints/round1-axis-local-inputs" \
  --matrix "$STEP/execution/checkpoints/round1-axis-local-formal" \
  --destination /tmp/arc-control-reliability-replay
```

No formal test expression is selected for scoring. Artifact verification hashes
the frozen source files without decoding test rows. The capped CPU run reproduced
all three CSVs byte for byte. Input and validation manifest digests are recorded in
[provenance](provenance.json).
The final run took 7.84 seconds and peaked at 2,434,352 KiB RSS under the
10 GiB address-space cap. It ran on CPU with two threads.
