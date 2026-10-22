# Real-model run cards (Stage 2 benchmark)

One directory per real adapter: `config.json` is the authoritative record of the
author defaults and the 18,164-gene setting; `README.md` is the run card
(environment, command, smoke result, and the full five-fold run table);
`FIDELITY.md` is the ticket-09 deviation ledger. The adapter code lives in
`../../adapters/<model>/run.py`.

| model | adapter | env | five-fold status | AUPRC@200 | ideal_fit |
|---|---|---|---|---|---|
| PRnet | `adapters/prnet/run.py` | `arc-stage2-bench-prnet` | 5/5 completed | 0.1643 | full |
| CPA | `adapters/cpa/run.py` | `arc-stage2-bench-cpa` | 5/5 completed | 0.0910 | full |
| chemCPA | `adapters/chemcpa/run.py` | `arc-stage2-bench-chemcpa` | 5/5 completed | 0.0608 | partial |
| biolord | `adapters/biolord/run.py` | `envs/biolord` | 5/5 completed | 0.0974 | partial |
| PerturbNet | `adapters/perturbnet/run.py` | `arc-stage2-bench-perturbnet` | 5/5 completed | **0.1861** | partial |
| Squidiff | `adapters/squidiff/run.py` | `arc-stage2-bench-squidiff` | 5/5 completed | 0.0212 | full (convergence failure) |
| CMonge | `adapters/cmonge/run.py` | `envs/cmonge` | 5/5 completed | 0.0732 | full |

The full 12-row leaderboard (baselines included), 95% CIs, secondary metrics,
the Squidiff convergence failure, and the recommendation are in
[`../REPORT.md`](../REPORT.md). No entrant beats the `delta_mean` baseline
(0.1989); PerturbNet is the best entrant.
