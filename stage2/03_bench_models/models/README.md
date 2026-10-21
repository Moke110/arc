# Real-model run cards (Stage 2 benchmark)

One directory per real adapter: `config.json` is the authoritative record of the
author defaults and the 18,164-gene setting; `README.md` is the run card
(environment, command, smoke result); `FIDELITY.md` is the ticket-09 deviation
ledger. The adapter code lives in `../../adapters/<model>/run.py`.

| model | adapter | env | status | ideal_fit |
|---|---|---|---|---|
| PRnet | `adapters/prnet/run.py` | `arc-stage2-bench-prnet` | fold-0 smoke complete | full |
| CPA | `adapters/cpa/run.py` | `arc-stage2-bench-cpa` | fold-0 smoke complete | full |
| chemCPA | `adapters/chemcpa/run.py` | `arc-stage2-bench-chemcpa` | fold-0 smoke complete | partial |
| biolord | `adapters/biolord/run.py` | `envs/biolord` | fold-0 smoke complete | partial |
| PerturbNet | `adapters/perturbnet/run.py` | `arc-stage2-bench-perturbnet` | fold-0 smoke complete | partial |
| Squidiff | `adapters/squidiff/run.py` | `arc-stage2-bench-squidiff` | fold-0 smoke complete | full |
| CMonge | `adapters/cmonge/run.py` | `envs/cmonge` | fold-0 smoke complete | full |

Ticket 16 runs all five cold-drug folds for these seven models plus the five CPU
baselines. Smoke numbers are wiring proofs, not rankings.
