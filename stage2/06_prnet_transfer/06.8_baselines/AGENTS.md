# AGENTS.md — 06.8_baselines (ticket #8)

Routing index (ADR 0030: ≤30 lines). Definitions, record and measurements live in
[`README.md`](README.md); the decision is [ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md)
and the frozen input is [`06.3_input`](../06.3_input/README.md).

- `execution/baselines.py` — the four baselines only: `no_change`, `delta_mean`
  (primary), `delta_mean_dose`, `ridge_fp` (auxiliary). Fits use **only** the fold's
  `role=train` treated responses; predictions are written in the `ptlib.predict`
  format and scored by `ptlib.evaluate.run_evaluation` — the one evaluator, never
  a re-implementation.
- `execution/tests/` — stdlib `unittest`, offline on the ticket-#4 fixture:
  fitters, leakage (poisoned valid/test targets), evaluator sanity on known
  ranking cases, record/determinism, and the bounded-memory block paths.
- `execution/tools/replicate_fixture.py` — test aid: row-replicated fixture for the
  memory gate (`make mem-check`); never evidence.
- `execution/Makefile`, `execution/run.sh` — `test`, `check`, `run`, `record`,
  `fixture`, `dry-run`, `mem-check`.
- `output/` — gitignored compute output with a narrow `.gitignore` exception for the
  tracked record only: `summary_table.csv`, `metrics.csv`, `per_drug_dose.csv`,
  `fold_coverage.csv`, `manifest.json`. Per-fold prediction/eval trees and fit
  params stay local.

Rules: never import or edit `06.4_pipeline/ptlib` internals to change a metric or
the truth Δ; never mutate the frozen 06.3 input; no `n_jobs`/worker pool; no
`observation × 18,164` dense array (gene-block reads, row-block writes,
`float32`); the real `--all` runs under `systemd-run … -p MemoryMax=10G` and
records its peak RSS in the ticket.
