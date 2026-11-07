# Tahoe pilot outcome

All six cadence profiles reached 1,000 optimizer steps. Their full
validation loaders evaluated 5,699, 5,698 and 5,700 groups in folds 0–2,
respectively. Each pair matched the author's EMA monitor exactly at steps
500 and 1,000. The denser cadence added two observations and achieved the
same strict-best loss in each fold.

| Fold | Cadence in steps | Validation seconds | Checkpoint seconds | Overhead relative to optimizer time | Best validation loss |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 250 | 42.70 | 11.74 | 95.34% | 4.530767 |
| 0 | 500 | 21.03 | 5.99 | 47.71% | 4.530767 |
| 1 | 250 | 42.35 | 12.56 | 97.83% | 4.518877 |
| 1 | 500 | 20.81 | 6.06 | 47.97% | 4.518877 |
| 2 | 250 | 41.71 | 11.86 | 94.83% | 4.421935 |
| 2 | 500 | 20.50 | 7.32 | 50.76% | 4.421935 |

The frozen rule selected 500 in every fold. The 250-step overhead exceeded
the protocol's 30% ceiling, and the 500-step overhead still exceeded the
15% target. Peak worker RSS across the profiles ranged from 4,574,056 to
5,097,628 KiB. Peak PyTorch allocation was 3,901,028,352 bytes per worker,
with 4,104,126,464 bytes reserved. The scheduler sampled up to 4.861 GiB
host process-tree RSS and 6.272 GiB of workstation GPU memory in use.

## Early stopping

All three calibrations reached 8,500 steps and reproduced their
500-step profile prefix. The candidate replay selected a 4,000-step
minimum in all folds. Each replay ended before the calibration horizon
and kept the best loss within 0.25% of the full trace's best.

| Fold | Patience in steps | Minimum improvement | Simulated stop | Best loss through stop | Full-trace best loss and step |
| --- | ---: | ---: | ---: | ---: | --- |
| 0 | 2,000 | 0.013968210 | 8,000 | 0.721680 | 0.720065 at 8,500 |
| 1 | 1,000 | 0.013911695 | 7,500 | 0.693151 | 0.693151 at 7,500 |
| 2 | 1,000 | 0.003445470 | 8,000 | 0.668565 | 0.667838 at 8,500 |

Calibration fit wall times were 714–733 seconds, with peak worker RSS
between 6,103,368 and 6,190,964 KiB. The scheduler sampled up to 5.904 GiB
host process-tree RSS and 6.204 GiB of workstation GPU memory in use.
The replayed stop steps describe the seed-17 calibration traces.
Formal fits reapply each fold's policy to its own validation losses.

[The selection manifest](selection/selection.json) freezes the three
fold policies and nine seed-specific formal recipes. The [protocol](PROTOCOL.md)
provides the replay commands. The [cadence evidence](cadence.json) retains
the validation histories, measured costs and exact input and recipe identities.
The [dispatch record](dispatch_identity.json) identifies the profile planner,
trainer, protocol and release at launch. The [calibration plan](calibration/plan.json)
binds the trainer and the paired traces to calibration. The [execution record](execution.json)
identifies the raw scheduler telemetry and completed attempts under
`execution/output/tahoe-pilot-ticket71-rng-fixed/`.
