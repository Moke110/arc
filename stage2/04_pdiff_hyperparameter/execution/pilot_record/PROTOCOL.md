# Step 04 validation-only pilot protocol

Freeze this protocol before reading pilot outcomes. Apply it independently
to each fold using seed 17, random initialization, the same frozen predictor
request, and the same training and validation settings from the shared
PerturbDiff trainer. No held-out responses participate. The 250- and
500-step fits both run for 2,000 steps with early stopping disabled. Each
uses the full validation loader at every scheduled event. The trainer
captures the native EMA monitor and performs its normal strict-best
checkpoint refreshes and reload. Pilot fits and replays are outside the formal
matrix.

## Cadence

Require both runs to complete their scheduled steps and pass the trainer's
reload check. Check the full validation batching contract and the monitor on
the shared 500-step grid. Reject a pair whose same-step losses differ beyond
the numeric tolerance in `protocol.json`. This checks that changing the
validation cadence did not change the underlying training trajectory.

For each cadence, sum its event validation and checkpoint seconds
separately, then divide their sum by optimizer training seconds. Record
setup, reload and wall time separately. Compare minimum observed monitor
values over the entire common step budget. Select 250 if its overhead
fraction is at most 15% of training and its best loss is within 0.25% of the
500-step best. Also select 250 if it improves the best loss by more than
0.25% and its overhead fraction is at most 30%. In all other cases choose
500 and record its measured overhead, including any over-target cost. The
15% target bounds routine checkpoint and validation cost. The 30% ceiling
allows a measurable gain in selection quality. The loss tolerance uses the
absolute 500-step best loss as its scale.

## Early stopping

The paired runs cover early optimization. For each fold, rerun at its chosen
cadence through 10,000 steps with stopping disabled. This horizon reaches into
the 6,000–13,000-step range of the previously reported Step 03 best checkpoints.
Verify that the calibration trace matches the chosen cadence's paired profile
at shared steps. The per-fit optimizer step and wall time caps and scheduler
grace are in `protocol.json`. A failed, incomplete, or mismatched trace cannot
freeze a recipe.

Replay all combinations of the three minimum training steps, three
optimizer-step patience windows, and three improvement thresholds in the
protocol through the trainer's stopping policy. Convert improvement
fractions to absolute minimum improvement using the magnitude of the fold's
first calibration monitor value. The candidate's score is its strict-best
native validation loss through its simulated stop. A candidate is eligible
only if it stops before the final calibration event and its score is within
0.25% of the full trace's strict-best validation loss. Pick the earliest
eligible stop, then the lowest score, then the smallest minimum steps,
smallest patience, and smallest improvement threshold, in that order.
This permits early stopping only when the observed optimizer savings preserve
the fold's validation result.

If no candidate is eligible, use the predeclared fallback of 10,000 minimum
steps, 4,000 patience steps and the smallest threshold. Mark the calibration
as right-censored and record whether policies reached the horizon or stopped
with excessive loss. The horizon and a single seed limit the strength of
convergence conclusions. Formal runs use the selected policy with the
fixed 30,000-step cap and seeds 17, 42 and 2026.

## Artifact workflow

Run the planner from the repository root. The output directory must be fresh.

```bash
pilot=stage2/04_pdiff_hyperparameter/execution/script/pilot.py
pilot_run=stage2/04_pdiff_hyperparameter/execution/output/tenv-pilot-reproduction
python "$pilot" profile "$pilot_run" --tenv-release stage2/03_modes_tenv/output/input
```

The dedicated scheduler runs `profile/queue.json`. After all six fits
complete, `python "$pilot" calibrate "$pilot_run"` validates the pairs,
writes `cadence.json` and freezes the calibration queue. The scheduler
then runs `calibration/queue.json`. Run `python "$pilot" select "$pilot_run"`
after the calibration queue completes. This emits `selection/selection.json`
and the nine trainer recipes beneath it. The manifest binds each fold input
identity and the effective recipe hash for every seed.

The compact evidence consists of `cadence.json`, `calibration/plan.json`,
and `selection/`. It preserves both profiles, the calibration monitor and
timings, the policy replays, scheduler peaks and raw trace hashes. Run
`python "$pilot" verify "$pilot_run"` with just that compact evidence to
replay both decisions and compare the emitted recipes.

## Probe random stream

The initial tenv pair diverged at matched optimizer steps. The upstream
DDIM prediction probe consumes global CUDA randomness even when its noise
coefficient is zero. Isolating probe draws from optimizer randomness
preserved the training stream across checkpoint refreshes. In the bounded
fold-0 regression, the 5- and 10-step cadence runs had identical validation
and training losses on the shared grid and identical selected EMA state
hashes. The initial pair is excluded from tuning and the bounded regression
is outside the pilot profile. See [probe evidence](tenv-probe-parity.json).
