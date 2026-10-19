# AGENTS.md — stage2

Routing index for Stage 2 (model setup and selection). Stage 1 lives at the
repository root in `stage1/` and is not a dependency here.

- `PLAN.md` — the plan: goal, setup stage, environments, 978 quarantine, exit
  criteria, deferred selection decisions.
- `01_setup/` — one run card per deployed model (`prnet`, `cpa`, `chemcpa`,
  `biolord`, `perturbnet`, `scvidr`, `squidiff`, `cmonge`, `decant`), each with
  its pinned environment, `run.sh`, and `output/`; `SUMMARY.md` records results.
- `audit/` — the 14-model entrant audit moved from `cm_bench_v1` (selection
  record; not an active benchmark).
- `_archive/cm_bench_v1/` — archived old benchmark; not a dependency.

Current phase: **setup only**. No benchmark metric, fold, or CM adapter is
defined here yet. Preprint models are excluded and their code is not vendored.
