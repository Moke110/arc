# Dose crossing with fixed tenv assessment parents

[Issue #264](https://github.com/Moke110/arc/issues/264) · November 21, 2026

The study crossed held doses 0.2, 1 and 10 µM within Step 03 cold-drug
outer fold 0 development rows. The outer held parents were excluded at
all doses. Within the frozen validation group, parent IDs were sorted by
SHA-256 and alternated between selection and assessment, using the same
rules as [#263](../optimization-selection.md). Each arm's assessment
used the same 15 parents, with 47 treated wells at the held dose.

| Held dose (µM) | Training and selection doses (µM) | Assessment setting |
|---:|---|---|
| 0.2 | 1 and 10 | Lower endpoint extrapolation |
| 1 | 0.2 and 10 | Interpolation |
| 10 | 0.2 and 1 | Upper endpoint extrapolation |

Each arm used 1,674 fit profiles from 275 parents and 98 selection
profiles from 16 parents. Both roles were balanced per parent across
the three dose cells before withholding one dose. The realized freeze
already balanced the replicate counts, so this rule removed zero profiles.
The assessment used all wells from complete validation parents on each
held dose. Each arm excluded 94 wells from the assessment parents at the
two other doses and 886 wells from fit and selection parents at the held
dose. [Run evidence](run.json) includes exact row identity digests, counts
by dose, and per-fit selection histories.

Both initializations ran with native Adam at 1e−4, seeds 42 and 17, a
batch size of 64 and 20 epochs (540 updates). The selection criterion was
maximum parent-equal Δ PCC with the earliest epoch as the tie-breaker.
Each arm used its train-only log-dose transform. The full-transfer arm
loaded the same frozen Tahoe winner in all experiments and verified the
dose-embedding affine bridge before optimization.

Each score uses predicted and observed expression change from the
same matched control. The direction error is the fraction of genes
with differing Δ signs, treating zero as its own direction. The amplitude
error is the mean absolute difference in predicted versus observed Δ
magnitudes. RMSE additionally measures signed Δ error. Transfer−random
is computed on the intersection of finite profile pairs per metric.
Summaries average profiles within parent, parents with equal weight,
and the two fixed seeds with equal weight. The percentile intervals
resample matched parent clusters across the three held doses 4,000 times
(seed 264). These intervals are conditional on the fits and paired
parents, with seed-specific point estimates reported separately.

## Results

| Held dose (µM) | Transfer − random Δ PCC | 95% interval | Seed 42 | Seed 17 | Δ RMSE change | Direction error change | Amplitude error change |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | +0.01830 | [−0.02457, +0.06651] | +0.01037 | +0.02623 | +0.00035 | +0.03418 | +0.00568 |
| 1 | −0.00384 | [−0.04470, +0.04688] | −0.01560 | +0.00792 | +0.00869 | +0.03423 | +0.00685 |
| 10 | −0.00648 | [−0.04096, +0.03245] | −0.02729 | +0.01433 | +0.01148 | +0.03703 | +0.00504 |

All 15 matched parents and 47 profiles per dose and seed had defined
paired PCC. The 0.2 µM point estimate was positive for both seeds.
At 1 and 10 µM the seeds disagreed on the sign of the gain. The paired
PCC gain contrast versus 1 µM was +0.02214 [−0.03333, +0.07423] at
0.2 µM and −0.00264 [−0.03826, +0.02547] at 10 µM. This prototype
does not establish a stable dose-dependent reversal. Holding parents
fixed removes the change of parent composition, while changing the held
dose still changes both the assessment condition and the fit dose support.

Direction error increased by 3.42–3.70 percentage points, with intervals
above zero at all three doses. Δ RMSE increased at 10 µM with a
conditional interval above zero. Positive error changes indicate worse
predictions. Amplitude error increased at all three doses, with intervals
above zero at 0.2 and 1 µM. Most transfer fits selected the epoch cap,
so these results describe short-budget transfer rather than convergence.

[Parent × dose × seed pairs](parent_pairs.csv) · [Paired estimates and intervals](summary.csv)

## Reproduction

The command uses the accepted Step 03 input pair and Step 04 Tahoe
checkpoint. The run starts from a new directory and does not train on
or score any held rows of the frozen outer fold. For a metadata-only check,
provide `--inputs --output --design-only`.

```bash
STEP=stage2/04_tsg_transfer
PAIR=stage2/03_tsg_dose/execution/checkpoints/round1-axis-local-inputs
"$STEP/execution/capped_phase.sh" /tmp/arc-264-evidence dose-crossing \
  python "$STEP/execution/dose_crossing.py" --inputs "$PAIR/tenv" \
  --tahoe-inputs "$PAIR/tahoe_bulk" \
  --experiment "$STEP/execution/checkpoints/accepted/experiment" \
  --source pert_model/TranSiGen/source --output /tmp/arc-264-replay
```

The run took 1:44.08 at 2,167,368 KiB peak RSS with the 10 GiB
systemd cap and shared GPU lease. The scope recorded all 12 fits and
terminal status 0. A separate metadata replay reproduced every role
digest and observed count. All 90 parent-dose-seed rows reconciled
with the per-fit assessment and reported paired means.
