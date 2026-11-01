# AGENTS.md — stage2/04_tahoe_data

Routing index. Purpose, corpus numbers and reproduce commands live in
`README.md`; the contract is ADR 0029. This step is **post-download processing
only**: it downloads nothing and imports nothing.

- `execution/build_pseudobulk.py` — the only entry point (scan → aggregate →
  assemble → qa); `--check` prints the resolved root, input and output. It finds
  the root by walking up to `.git`; never add `parents[N]`/`../..` arithmetic.
- Input: the imported atlas `pert_data/Tahoe-100M/source/`; the importer and its
  verifiers live in `pert_data/Tahoe-100M/script/`.
- `output/pseudobulk/tahoe_pseudobulk.h5ad` — 25.8 GB corpus, local-only, gitignored.
- `output/pseudobulk/{qa.json,qa_report.md,qa_contrasts.csv,manifest.csv}` —
  tracked build record: shape, QA gates, per-contrast agreement, profile table.
- `output/pseudobulk/_parts/` — 14 per-plate CSR parts, the resumable part cache
  (intermediate, gitignored; `aggregate` rebuilds a missing part).
- `output/pseudobulk/_state/` — `shard_plate.json` scan bookkeeping (machine
  state, gitignored; safe to delete, `scan` rebuilds it).
- `--plates plate4,plate5` limits any stage; `--force` rebuilds existing parts.
- `execution/build_pretrain_view.py` — **planned, nothing built**: builds
  `output/pretrain_view/tahoe_cm_schema.h5ad` (67,018 x 18,164 CM-schema view)
  from `output/pseudobulk/_parts/`. Contract + deviations D0–D3:
  `.scratch/stage2-05-prnet-pretrain/spec.md`. The 2026-11-01 build reached
  assembly but never ran QA, and its unverified output was deleted the same day;
  `output/pretrain_view/` does not currently exist. Do not rebuild the corpus.

Reproduce commands: `README.md` (`--check`, `qa`, and the full rebuild).
