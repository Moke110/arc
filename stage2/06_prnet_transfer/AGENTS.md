# AGENTS.md — stage2/06_prnet_transfer

Routing index. Purpose, first-round contract, status and reproduce commands live in
[`README.md`](README.md); the design decision is
[ADR 0032](../../docs/adr/0032-prnet-formal-transfer.md). Frozen evidence in
`stage2/02_cm_data/output/`, `stage2/04_tahoe_data/output/` and
`stage2/05_prnet_pretrain/output/` is read-only here; every `output/` below is
gitignored compute output with a tracked record beside it.

- `eligibility/` — the 460-structure pharmacological eligibility audit (evidence
  table, ATC/MoA normalisation, coverage counts).
- `qc/` — per-observation experimental QC, control matching, thresholds locked
  before any model result, and the final observation/well counts.
- `input/` — frozen two-source 18,164-gene input: raw counts + `log1p_cp10k`,
  missing-gene mask, five drug-key folds, per-fold consensus-control membership,
  locked executable config.
- `pipeline/` — the only training/evaluation code: data loading, config-driven
  structure, two-stage Gaussian-NLL training with direction-AUPRC@200 selection,
  prediction, independent scoring, run manifests.
- `runs/` — three seeded Tahoe pretrains; 45 random-init and 45 Tahoe-init CM
  candidates; selection records.
- `baselines/` — the five-fold constant / training-mean baselines on the identical
  test population and evaluator.
- `evaluation/` — the 30 selected checkpoints, paired gains, stratification, model
  cards and the report.
- `ablations/capacity/`, `ablations/env_compounds/` — post-first-round ablations,
  explicitly not first-round completion gates.

Rules: find the repo root by walking up to `.git` (never `Path.parents`); never
rewrite frozen evidence or a `source/` tree; `pipeline/` is the only place
training/evaluation code lives; selection never reads test treated expression.
