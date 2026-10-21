# AGENTS.md — stage2

Routing index for Stage 2 (model setup and selection). Stage 1 lives at the
repository root in `stage1/` and is not a dependency here.

- `PLAN.md` — the plan: goal, setup stage, environments, 978 quarantine, exit
  criteria, deferred selection decisions.
- `01_setup/` — one run card per deployed model (`prnet`, `cpa`, `chemcpa`,
  `biolord`, `perturbnet`, `scvidr`, `squidiff`, `cmonge`, `decant`), each with
  its pinned environment, `run.sh`, and `output/`; `SUMMARY.md` records results.
- `02_cm_data/` — Stage 2 CM perturbation corpus and dataset bundles: frozen
  `corpus/v1/` (18,164 genes, 525 drug keys) and built `datasets-1/` (five
  h5ad views). Build contract/results in `README.md`; decisions in ADR 0026.
- `03_bench_models/` — executed benchmark: frozen ticket-11 `data/`, the
  `adapter_api`/`runner`/`evaluate` harness, five CPU baselines, and all seven
  real adapters (`prnet`, `cpa`, `chemcpa`, `biolord`, `perturbnet`,
  `squidiff`, `cmonge`) with configs/run cards/ledgers in `models/`; ticket 16
  runs the five-fold queue. See `README.md`.
- `audit/` — the 14-model entrant audit moved from `cm_bench_v1` (selection
  record; not an active benchmark).
- `_archive/cm_bench_v1/` — archived old benchmark; not a dependency.

Current phase: **benchmark build**. `03_bench_models/` has the frozen data,
harness, CPU baselines, and all seven real adapters (each fold-0 smoke-tested);
the five-fold run is the remaining execution work. Preprint models are excluded
and their code is not vendored.
