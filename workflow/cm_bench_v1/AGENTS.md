# AGENTS.md — workflow/cm_bench_v1

Routing index for the Stage 2 v1 whole-transcriptome cardiomyocyte benchmark
(`cm_bench_v1`). Stage 2 decides which open-source whole-transcriptome
drug-perturbation model should produce the perturbomes that feed future MIRI
scoring (ADR 0020); it produces no MIRI drug ranking. One line per entry.

- `run.sh` — end-to-end: audit → `01_bench_data` → `04_evaluate`.
- `audit_entrants.py` — pre-flight checkpoint-availability gate (ticket 02).
- `01_bench_data/` — raw CM datasets → per-plate Δ → common gene universe → limma
  true sets → cold-drug/cold-dose masks; gate assertions (tickets 01/03/04).
- `02_baselines/baselines.py` — the four mandatory baselines implementing the
  `predict` seam (ticket 05); `02_baselines/gate.py` is the in-sample sanity gate.
- `03_<model>/` — one adapter directory per model entrant (tickets 07-11);
  `03_stub/` proves the runner harness (entry + adapter + gate).
- `04_evaluate/evaluate.py` — direction-aware AUPRC, aggregation, leaderboard,
  report (tickets 01/04/12).
- `_lib/` — shared: `adapter.py` (the one seam), `model_runner.py` (per-env
  subprocess bridge), `model_registry.py`, `metric.py`, `true_sets.py`,
  `splits.py`, `genes.py`, `datasets.py`, `normalize.py`, `baseline.py`,
  `bench_io.py`, `gates.py`, `evaluate.py`, `entrant_audit.py`,
  `primary_concentrations.py`, `utils.py`, and the `fetch_*.py` asset builders.
- `input/` — frozen assets: `ensembl_to_symbol.tsv`, `gse244740_smiles.tsv`.
- `output/` — `bench/` artifacts, leaderboards, `report.md`, `entrant_audit.csv`.
- `docs/entrant_audit.md` — the trimmed, documented entrant list.

## The one seam

Every entrant implements `predict(train_cases, test_cases) -> Δ` (`_lib/adapter.py`).
`04_evaluate` and the metric are exercised end-to-end on CPU through the
`no-change` adapter, so the evaluation path is testable with no model env or GPU.
Models are reached only through this seam, via `_lib/model_runner.py` (their own
pinned conda env, ADR 0021) — never by importing model internals.

## Scope guard

No 978-landmark artifact is read anywhere (ADR 0020). `gse262419` and any in-HG
drug restriction are excluded. Stage 1 workflows/components/data are untouched.
