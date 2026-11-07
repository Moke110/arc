# Tahoe validation-only pilot

[protocol.json](protocol.json) fixes the Tahoe budgets and candidate grid.
The shared [pilot comparison procedure](../pilot_record/PROTOCOL.md) governs
cadence and stopping selection using these Tahoe parameters. Freeze the
protocol before dispatch and use it independently within each cold-drug
fold. Use seed 17 and keep the model, learning rate, scheduler, batch size
and input fixed. Only training and validation expression enter the pilot.
The held-out covariate vocabulary is part of the frozen predictor request.

## Budget and comparison

The two random-init profiles per fold each run to 1,000 optimizer steps.
This gives four full validation events at cadence 250, two at cadence 500,
and two shared steps for checking that cadence leaves the training trajectory
unchanged. The independent profile cap bounds the cost measurement. Selection
uses the native EMA monitor and normal checkpoint refreshes, retaining
validation time, checkpoint time, worker peak RSS, and GPU memory peaks.

Use the measured sum of validation and checkpoint seconds divided by
optimizer training seconds. Choose cadence 250 when it is within 0.25% of
the best loss at cadence 500 and its overhead is at most 15% of training.
A gain greater than 0.25% permits up to 30% overhead. Otherwise use cadence
500 and record its overhead. Both profiles must complete and agree at shared
steps within the numeric tolerance in the protocol.

At eight training groups per batch, the verified Tahoe release has 3,992,
3,992 and 4,012 optimizer steps per pass through the training groups in
folds 0, 1 and 2, respectively. The 8,500-step calibration includes two
full passes in each fold. This horizon gives the latest candidate minimum
of 6,000 steps room for its full 2,000-step patience window and one further
500-step validation event. Earlier minimums near and after the first pass
are 4,000 and 5,000 steps. The single-seed trace screens for early policies
that preserve the observed loss across the calibration horizon. Use the
predeclared right-censored fallback if every candidate reaches the horizon
or loses too much validation performance. The formal step cap is 30,000.

Calibrate using a fresh run at the selected cadence and verify agreement
with the short profile over their shared steps. Replay the shared trainer's
step-based policy over the entire candidate grid using the first monitor
value's absolute magnitude to scale each improvement threshold. Apply the
shared loss tolerance and choose the earliest eligible stop, then lower
loss, then smaller minimum steps, patience, and threshold. Each event can
refresh the strict-best checkpoint regardless of the stopping threshold.
The pilot output is separate from the formal matrix.

## Reproduce the queue

From the repository root, use the verified Tahoe release from the full-axis
Tahoe smoke. Use a fresh output path and give each generated queue to the
dedicated scheduler.

```bash
step=stage2/04_pdiff_hyperparameter
execution="$step/execution"
pilot="$execution/script/pilot.py"
record="$execution/output/tahoe-pilot-replay"
release="$execution/output/tahoe-smoke-ticket68/release"
python "$pilot" profile "$record" --tahoe-release "$release" \
    --protocol "$execution/tahoe_pilot_record/protocol.json"
# Run "$record/profile/queue.json" through the dedicated scheduler.
python "$pilot" calibrate "$record"
# Run "$record/calibration/queue.json" through the dedicated scheduler.
python "$pilot" select "$record"
python "$pilot" verify "$record"
```

The compact cadence evidence is in `cadence.json`. The verified selection
manifest and its nine formal recipe files are in `selection/`.
