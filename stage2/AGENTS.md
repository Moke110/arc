# AGENTS.md — stage2

Routing index for Stage 2. One line per step; purpose, results, constraints and
reproduce commands live in each step's `README.md`; overview: [`README.md`](README.md).
Stage 1 is `stage1/`, not a dependency here; the data contract is ADR 0030.

- [`01_setup/`](01_setup/README.md) — nine deployed non-preprint models proven to
  load/train in pinned per-model environments; run cards under `execution/<model>/`.
- [`02_dataset_corpus/`](02_dataset_corpus/README.md) — every source dataset as a
  native-axis raw-count h5ad corpus: five corpora in `output/corpus/v2/`, seven
  crosswalk tables + coverage README, `--verify` G1–G6 PASS (ADR 0033, amends 0026).
- [`_archive/03_models_cmbench/`](_archive/03_models_cmbench/README.md) — executed 60-unit cold-drug
  benchmark, **frozen pending a redo on the merged corpus**; inputs now point at
  `stage2/_frozen/cm-data-v1/` (ADR 0027, ADR 0028).
- [`03_modes_tenv/`](03_modes_tenv/README.md) — active shared-Morgan tenv comparison;
  frozen input, discovery, runtime admission and artifact pipeline (ADR 0042).
- **`04` is a deliberate gap** — the retired Tahoe data step was folded into
  `02_dataset_corpus/` (ADR 0033); the number is not reused.
- [`05_prnet_pretrain/`](05_prnet_pretrain/README.md) — executed PRnet pretrain
  (2026-11-01; pipeline claim only) and owner of the CM-schema Tahoe view + builder
  it trains on (ADR 0031, ADR 0033).
- [`06_prnet_transfer/`](06_prnet_transfer/README.md) — formal-transfer study
  (executed 2026-11-02, paths only reconnected after the merge): frozen five-fold
  input, 3 Tahoe pretrains × 90 CM candidates, 30-unit paired evaluation; status
  **negative** first round, ablations #10/#11 deferred (ADR 0032).
- [`_envs/`](_envs/README.md) — non-step: shared per-model environments + locks.
- [`audit/`](audit/README.md) — non-step: entrant/fairness audit, 2026-10-26
  selection decisions.
- [`_archive/`](_archive/README.md) — non-step: superseded `cm_bench_v1/` and the
  archived `PLAN.md`; nothing here is an active dependency.
