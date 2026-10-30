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
  `03_bench_models/temp_env_all_5k/` is a paired sensitivity fork with the same
  wells, splits, entrants and seeds on a 5,000-gene HVG axis (rank-sum over
  `dtoxs_cm`, `screen_ctox`, `temp_env_all-smiles`); decision and leakage
  declaration in ADR 0028. It is not a replacement leaderboard.
- `04_tahoe_pretrain/` — Tahoe-100M data home + the Stage 2 bulk pretrain corpus: the full atlas (3,388 shards, 95.6M cells, 428.8 GB) is imported under `pert_data/Tahoe-100M/source/`; `script/build_pseudobulk.py` builds the well-level `(plate, well, cell_line)` pseudobulk corpus (67,018 profiles × 62,710 genes, raw counts + `log1p_cp10k`, `DMSO_TF` as control rows) into `output/pseudobulk/`, with `qa_report.md` vs the authors' DESeq2 DE; contract + measurements in ADR 0029; the importer lives in the item (`pert_data/Tahoe-100M/script/{fetch_tahoe.py,verify_import.py}`).
- `audit/` — the 14-model entrant audit moved from `cm_bench_v1` (selection
  record; not an active benchmark), the 2026-10-26 fairness audit, and
  `model_selection_decisions.md` (governs current follow-up candidates).
- `_archive/cm_bench_v1/` — archived old benchmark; not a dependency.

Current phase: **model selection evidence**. `03_bench_models/temp_env_all_18k/`
has the frozen data, harness, CPU baselines, and all seven real adapters,
executed over five cold-drug folds (60/60 units).
`temp_env_all_18k/REPORT.md` records that no entrant beats the `delta_mean`
baseline on the primary metric; PerturbNet is the best entrant.
`temp_env_all_5k/` repeats the run on a 5,000-gene HVG axis to test whether that
ranking depends on the transcriptome width (ADR 0028): PRnet moves above the
`delta_mean` baseline there, while the 18k best entrant PerturbNet falls behind,
so the 18k "no entrant beats the baseline" conclusion is gene-space dependent.
Preprint models are excluded and their code is not vendored.

**Current selection override (2026-10-26, user-approved):** Squidiff is excluded
from future selection, tuning, and debugging investment after both the 18k and
5k recipes failed to converge; further remediation has uncertain benefit.
Remaining candidates are PRnet, CPA, chemCPA, biolord, PerturbNet, and CMonge,
with a 15,000+ gene delivery requirement. Retain all historical seven-model
artifacts and registries for reproduction. New queues exclude Squidiff unless
the user reopens the decision; see `audit/model_selection_decisions.md`.
