# tenv pilot outcome

All six cadence profiles reached 2,000 optimizer steps. The two runs for
each fold produced identical author-monitor losses at their four shared
validation steps. The 500-step cadence kept the same best validation loss
while reducing validation and checkpoint time relative to optimizer time.

| Fold | Best validation loss at either cadence | 250-step overhead | 500-step overhead |
| --- | ---: | ---: | ---: |
| 0 | 1.242252 | 22.98% | 11.34% |
| 1 | 1.291144 | 21.75% | 10.81% |
| 2 | 1.260600 | 21.56% | 10.96% |

Each selected-cadence calibration reached 10,000 optimizer steps and
reproduced its 2,000-step profile prefix. All folds selected a minimum of
4,000 steps and patience of 1,000 steps. The absolute threshold reflects
each fold's initial validation loss and its replayed candidates.

| Fold | Minimum improvement | Simulated stop | Best loss through stop | Full-trace best loss and step |
| --- | ---: | ---: | ---: | --- |
| 0 | 0.005826172 | 6,500 | 0.868801 | 0.867684 at 7,500 |
| 1 | 0.005903126 | 5,500 | 0.866392 | 0.864349 at 7,000 |
| 2 | 0.002887517 | 7,000 | 0.832903 | 0.832903 at 6,500 |

The three calibration workers took 618–628 seconds each, with peak RSS
between 3,183,472 and 3,184,008 KiB. The cadence and calibration evidence
records separate validation and checkpoint timings, worker and workstation
peaks, resource caps, frozen input and trainer identities, and trace hashes.

The [selection manifest](selection/selection.json) freezes the three fold
policies and nine seed-specific formal recipes. The calibration used seed
17 and a 10,000-step horizon. Its replayed stop steps describe those
pilot traces. Formal fits reapply the frozen policies to their own traces.
The [protocol](../PROTOCOL.md) provides the replay command. Raw evidence
is local under `execution/output/tenv-pilot-ticket70-rng-isolated/`.
