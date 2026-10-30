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

Reproduce commands: `README.md` (`--check`, `qa`, and the full rebuild).
