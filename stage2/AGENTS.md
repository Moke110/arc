# AGENTS.md — stage2

Routing index for Stage 2 (model setup and selection). One line per step; purpose,
results, selection constraints and reproduce commands live in each step's
`README.md`, and the stage overview is [`README.md`](README.md). Stage 1 is
`stage1/`, not a dependency here; the data contract is ADR 0030.

- [`01_setup/`](01_setup/README.md) — nine deployed non-preprint models proven to
  load/train in their own pinned environments; run cards under `execution/<model>/`
  (ADR 0021, ADR 0024).
- [`02_cm_data/`](02_cm_data/README.md) — frozen model-agnostic CM corpus
  (`output/corpus/v1/`) and the five `output/datasets-1/` h5ad views (ADR 0026).
- [`03_models_cmbench/`](03_models_cmbench/README.md) — the executed 60-unit
  cold-drug benchmark in two substeps: `03.1_cold_drug_18k/` (primary, ADR 0027)
  and `03.2_hvg_sensitivity_5k/` (HVG sensitivity fork, ADR 0028).
- [`04_tahoe_data/`](04_tahoe_data/README.md) — the Tahoe-100M pseudobulk
  pretrain corpus built from the already-downloaded atlas (ADR 0029), plus its
  CM-schema pretrain view: 67,018 x 18,164, built and verified 2026-11-01, nine
  gates 9/9 PASS (ADR 0031).
- [`05_prnet_pretrain/`](05_prnet_pretrain/README.md) — the executed PRnet
  pretrain on 04's view (2026-11-01): 68 epochs, early stop at epoch 68, best
  valid Gaussian NLL −2.628080 at epoch 47, 883.2 s on GPU; artifacts under
  `output/` (ADR 0031). Pipeline claim only.
- [`_envs/`](_envs/README.md) — non-step: shared per-model benchmark
  environments and their `pip freeze` locks (ADR 0021).
- [`audit/`](audit/README.md) — non-step: entrant audit, fairness audit, and the
  2026-10-26 selection decisions.
- [`_archive/`](_archive/README.md) — non-step: superseded `cm_bench_v1/`
  benchmark and the archived `PLAN.md`; nothing here is an active dependency.
