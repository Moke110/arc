# AGENTS.md — stage2/06_prnet_transfer

Routing index for the PRnet formal-transfer step. Purpose, the first-round
contract, status and reproduce commands live in [`README.md`](README.md); the
research decisions live in the approved spec
(`.scratch/stage2-prnet-formal-transfer/spec.md`, gitignored) and the
architecture decision record written by the input ticket.

Everything under `output/` is gitignored compute output: a fresh clone
reproduces it from the committed `execution/` code, the tracked summary tables
and the manifests. Frozen evidence in `stage2/02_cm_data/output/`,
`stage2/04_tahoe_data/output/` and `stage2/05_prnet_pretrain/output/` is
**read-only** here.

- `eligibility/` — the 460-structure pharmacological eligibility audit: evidence
  table (include / exclude / uncertain), ATC/MoA normalisation, coverage counts.
- `qc/` — per-observation experimental QC audit for Tahoe and temp_env: raw
  plate/well provenance, control matching, quality flags locked before any model
  result, and the final observation/well counts.
- `input/` — the frozen two-source 18,164-gene input contract: raw counts +
  `log1p_cp10k`, missing-gene mask, the five drug-key folds, explicit per-fold
  consensus-control membership, and the locked executable config.
- `pipeline/` — the reusable library and CLIs: data loading, structure/recipe
  configs, two-stage training (Gaussian NLL), direction-AUPRC@200 checkpoint
  selection, prediction, independent scoring, and the run manifest writer.
- `runs/` — the formal executions: three seeded Tahoe pretrains, 45 random-init
  CM candidates, 45 Tahoe-init CM candidates, plus their selection records.
- `baselines/` — the five-fold constant / training-mean baselines evaluated on
  the identical test population with the identical evaluator.
- `evaluation/` — the 30 selected checkpoints, paired pretraining-gain and
  baseline comparisons, stratified analyses, model cards and the report.
- `ablations/capacity/`, `ablations/env_compounds/` — the post-first-round
  ablations (state latent 128 / hidden 256; extra environment compounds). These
  are explicitly **not** first-round completion gates.

Rules:

- `repo_root.py`-style root discovery (walk up to `.git`) is the only way to
  locate the repository; never use `Path.parents` or `../..` arithmetic.
- Never rewrite the frozen `stage2/02`/`04`/`05` outputs, and never edit a
  `source/` tree. New inputs, splits, configs and results get their own version
  identity.
- `pipeline/` is the only place training/evaluation code lives; `runs/`,
  `baselines/`, `evaluation/` and `ablations/` invoke it and record results.
- Selection never reads test-split treated expression; the evaluator is the only
  component that reads it.
