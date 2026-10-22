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
- `_envs/` — shared Stage 2 per-model bench environments (four
  `arc-stage2-bench-*` conda envs by name, `biolord`/`cmonge` venvs) and their
  `pip freeze` locks; promoted from the benchmark so later stages can reuse them.
- `03_bench_models/` — executed cold-drug benchmark, fully archived under
  `temp_env_all_18k/`: frozen `data/`, harness `scripts/` (original structure),
  per-unit `output/`, `results/`, configs, run cards, fidelity ledgers,
  `README.md`, `REPORT.md`, plus `summary.png` / `summary_table.csv` (the
  consolidated results table and per-metric bars). Five CPU baselines and all
  seven real adapters (`prnet`, `cpa`, `chemcpa`, `biolord`, `perturbnet`,
  `squidiff`, `cmonge`); 60/60 units completed. No entrant beats the
  `delta_mean` baseline; PerturbNet is the best entrant.
- `audit/` — the 14-model entrant audit moved from `cm_bench_v1` (selection
  record; not an active benchmark).
- `_archive/cm_bench_v1/` — archived old benchmark; not a dependency.

Current phase: **model selection evidence**. `03_bench_models/temp_env_all_18k/`
has the frozen data, harness, CPU baselines, and all seven real adapters,
executed over five cold-drug folds (60/60 units).
`temp_env_all_18k/REPORT.md` records that no entrant beats the `delta_mean`
baseline on the primary metric; PerturbNet is the best entrant.
Preprint models are excluded and their code is not vendored.
