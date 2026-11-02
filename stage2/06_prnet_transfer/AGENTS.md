# AGENTS.md — stage2/06_prnet_transfer

Routing index. Purpose, first-round contract, status and reproduce commands live in
[`README.md`](README.md); the design decision is
[ADR 0032](../../docs/adr/0032-prnet-formal-transfer.md). Layout follows
[ADR 0030](../../docs/adr/0030-external-data-provenance-and-step-layout.md):
`execution/<area>/` is harness/code, `output/<area>/` is the artifacts a later
ticket reads. Frozen evidence in `stage2/02_cm_data/output/`,
`stage2/04_tahoe_data/output/` and `stage2/05_prnet_pretrain/output/` is read-only
here; every `output/` is gitignored and carries a tracked record (`git add -f`).

- `execution/eligibility/` → `output/eligibility/` — the 460-structure
  pharmacological eligibility audit (evidence table, ATC/MoA normalisation,
  coverage counts, handoff set).
- `execution/qc/` → `output/qc/` — per-observation experimental QC, control
  matching, thresholds locked before any model result, final counts, handoff.
- `execution/input/` → `output/input/` — frozen two-source 18,164-gene input,
  five drug-key folds, per-fold consensus-control membership, locked config.
- `execution/pipeline/` — the only training/evaluation code: config-driven
  structure, two-stage Gaussian-NLL training with direction-AUPRC@200 selection,
  prediction, independent scoring, run manifests.
- `execution/runs/` → `output/runs/` — three seeded Tahoe pretrains; 45
  random-init and 45 Tahoe-init CM candidates; selection records.
- `execution/baselines/` → `output/baselines/` — five-fold constant /
  training-mean baselines on the identical test population and evaluator.
- `execution/evaluation/` → `output/evaluation/` — the 30 selected checkpoints,
  paired gains, stratification, model cards and the report.
- `execution/ablations/{capacity,env_compounds}/` → `output/ablations/…` —
  post-first-round ablations, explicitly not first-round completion gates.

Rules: find the repo root by walking up to `.git` (never `Path.parents`); never
rewrite frozen evidence or a `source/` tree; `execution/pipeline/` is the only
place training/evaluation code lives; selection never reads test treated
expression.
