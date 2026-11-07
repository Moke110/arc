# AGENTS.md — stage2

Routing index for Stage 2. One line per step; purpose, results, constraints and
reproduce commands live in each step's `README.md`; overview: [`README.md`](README.md).
Stage 1 is `stage1/`, not a dependency here; the data contract is ADR 0030.

- [`01_setup/`](01_setup/README.md) — native per-model setup cards and pinned
  environments, including PerturbDiff; run cards under `execution/<model>/`.
- [`02_dataset_corpus/`](02_dataset_corpus/README.md) — every source dataset as a
  native-axis raw-count h5ad corpus: five corpora in `output/corpus/v2/`, seven
  crosswalk tables + coverage README, `--verify` G1–G6 PASS (ADR 0033, amends 0026).
- [`03_modes_tenv/`](03_modes_tenv/README.md) — completed shared-Morgan tenv
  comparison on the 18,164-gene CM axis (ADR 0042).
- `04_pdiff_tahoe` and `05_pdiff_tenv` — planned standalone PerturbDiff
  training studies. See ADR 0046 for the agreed scope.
- [`_envs/`](_envs/README.md) — non-step: shared per-model environments + locks.
- [`audit/`](audit/README.md) — non-step: entrant/fairness audit, 2026-10-26
  selection decisions.
- [`_archive/`](_archive/README.md) — non-step: superseded `cm_bench_v1/` and the
  archived `PLAN.md`; nothing here is an active dependency.
