# AGENTS.md — stage2/04_tahoe_data

Routing index. Purpose, numbers and reproduce commands live in `README.md`; the
corpus contract is ADR 0029, the view's is
`.scratch/stage2-05-prnet-pretrain/spec.md`. Post-download processing only.

- `execution/build_pseudobulk.py` — corpus: scan → aggregate → assemble → qa;
  `--check` resolved paths; `--plates` limits a stage, `--force` rebuilds parts.
  Root by walking up to `.git`, never `parents[N]`/`../..`; input is the imported
  atlas `pert_data/Tahoe-100M/source/` (importer + verifiers in its `script/`).
- `output/pseudobulk/tahoe_pseudobulk.h5ad` — 25.8 GB corpus, local-only, ignored.
- `output/pseudobulk/{qa.json,qa_report.md,qa_contrasts.csv,manifest.csv}` — the
  tracked corpus build record (shape, gates, per-contrast agreement, profiles).
- `output/pseudobulk/_parts/` — 14 per-plate CSR parts, resumable cache
  (gitignored; `aggregate` rebuilds a missing part); `_state/shard_plate.json` is
  scan bookkeeping (safe to delete, `scan` rebuilds it).
- `execution/build_pretrain_view.py` — view, built + verified 2026-11-01 (9/9
  gates), writes `output/pretrain_view/`. Default = produce + QA; `--qa` =
  validate only (non-zero on failure); `--check` = paths; `--force` = re-gather.
  Writer `arc` (anndata 0.12.16), reader `arc-stage2-bench-prnet` (0.10.8).
- `output/pretrain_view/tahoe_cm_schema.h5ad` — 67,018 x 18,164 CSR int32 +
  `log1p_cp10k`, 16.32 GB, local-only, gitignored.
- `output/pretrain_view/{manifest.csv,compounds.csv,qa.json,qa_report.md}` — the
  tracked view record (counts + digest, 379 compounds, gates, report). PRnet was
  pretrained on the view 2026-11-01 (68 epochs, best valid NLL −2.628080; ADR 0031).
- `output/pretrain_view/_view_parts/` — per-plate gather cache (gitignored
  intermediate; a missing part is rebuilt, `--force` re-gathers all).
- The corpus builder is **not** a dependency of the view build: never imported
  or run, no shard re-aggregated — the view reads the frozen obs + `_parts/`.
