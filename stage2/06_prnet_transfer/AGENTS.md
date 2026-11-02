# AGENTS.md — stage2/06_prnet_transfer

Routing index. Study design, fixed contract and status live in
[`README.md`](README.md); the study decision is
[ADR 0032](../../docs/adr/0032-prnet-formal-transfer.md). Per ADR 0030 §3 this
step is a container: each of the eleven substeps below carries its own
`AGENTS.md`, `README.md`, `execution/` and `output/`, and every `output/` is
gitignored compute output with a tracked record force-added beside it. Frozen
evidence in `stage2/{02,04,05}_*/output/` and every `source/` tree is read-only.

- [`06.1_eligibility/`](06.1_eligibility/README.md) — #1: 460-structure audit (include/exclude/uncertain, ATC+MoA normalisation, coverage counts, handoff).
- [`06.2_qc/`](06.2_qc/README.md) — #2: per-observation QC, control matching, thresholds locked before any model result, counts and handoff.
- [`06.3_input/`](06.3_input/README.md) — #3: frozen two-source 18,164-gene input, five drug-key folds, per-fold consensus controls, locked config.
- [`06.4_pipeline/`](06.4_pipeline/README.md) — #4: the only training/evaluation code; real-data smoke of both initialisations, AUPRC@200 selection, scoring.
- [`06.5_tahoe_pretrain/`](06.5_tahoe_pretrain/README.md) — #5: three seeded formal pretrains, validation curves, model cards.
- [`06.6_cm_random/`](06.6_cm_random/README.md) — #6: 45 random-init CM candidates and the per-fold LR selection record.
- [`06.7_cm_tahoe/`](06.7_cm_tahoe/README.md) — #7: 45 Tahoe-init CM candidates with strict, axis-checked warm start.
- [`06.8_baselines/`](06.8_baselines/README.md) — #8: five-fold constant / training-mean baselines on the identical test population and evaluator.
- [`06.9_evaluation/`](06.9_evaluation/README.md) — #9: 30 selected checkpoints, paired pretraining gains, baseline comparison, stratification, report.
- `06.10_capacity_ablation/` (not built) — #10, deferred: state latent 128 / hidden 256, post-first-round, not a completion gate.
- `06.11_env_compound_ablation/` (not built) — #11, deferred: extra environment compounds, post-first-round, not a completion gate.

Rules: find the repo root by walking up to `.git` (never `Path.parents`); never
rewrite frozen evidence; `06.4_pipeline/execution/` is the only place
training/evaluation code lives; selection never reads test treated expression.
