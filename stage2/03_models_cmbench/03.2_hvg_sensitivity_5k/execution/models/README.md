# Real-model run cards (Stage 2 5k sensitivity benchmark)

One directory per real adapter: `config.json` is the authoritative record of the
author defaults and the 5,000-gene setting; `README.md` is the run card
(environment, command, smoke result, and the full five-fold run table);
`FIDELITY.md` is the ticket-09 deviation ledger, re-declared for the 5k axis.
The adapter code lives in `../scripts/adapters/<model>/run.py`.

This is the sensitivity fork of the executed 18k benchmark; architecture,
optimizer, schedule, seed and budget are unchanged. Only the output width and
the gene-space rows of each ledger differ.

| model | adapter | env | five-fold status | AUPRC@200 | gap vs `delta_mean` | ideal_fit |
|---|---|---|---|---|---|---|
| PRnet | `scripts/adapters/prnet/run.py` | `arc-stage2-bench-prnet` | 5/5 completed | **0.2735** | **+0.0243** | full |
| PerturbNet | `scripts/adapters/perturbnet/run.py` | `arc-stage2-bench-perturbnet` | 5/5 completed | 0.1929 | −0.0563 | partial |
| CPA | `scripts/adapters/cpa/run.py` | `arc-stage2-bench-cpa` | 5/5 completed | 0.1559 | −0.0932 | full |
| biolord | `scripts/adapters/biolord/run.py` | `stage2/_envs/biolord` | 5/5 completed | 0.1467 | −0.1025 | partial |
| chemCPA | `scripts/adapters/chemcpa/run.py` | `arc-stage2-bench-chemcpa` | 5/5 completed | 0.1165 | −0.1327 | partial |
| CMonge | `scripts/adapters/cmonge/run.py` | `stage2/_envs/cmonge` | 5/5 completed | 0.1153 | −0.1339 | full |
| Squidiff | `scripts/adapters/squidiff/run.py` | `arc-stage2-bench-squidiff` | 5/5 completed | 0.0523 | −0.1968 | full (convergence failure) |

`delta_mean` scores 0.2491 on the same folds. PRnet is the only entrant above
that baseline; its per-drug paired CI is +0.0113..+0.0367. The full 12-row
leaderboard, the paired 18k comparison and the Squidiff convergence failure are
in [`../../README.md`](../../README.md). The executed 18k result remains primary.
