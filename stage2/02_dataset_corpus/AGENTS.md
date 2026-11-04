# AGENTS.md — stage2/02_dataset_corpus

Routing index. Human record: `README.md`; contracts: ADR 0033 (the corpus
layer) and ADR 0026 (CM identity). Status: skeleton (ticket #14) — tickets
#15/#16 build the corpora/tables, #18 owns the downstream view, #20 docs.

- `execution/build_cm_datasets.py` — CLI today; `--check` prints the four CM
  source counts. Ticket #15 replaces it with a unified `build_corpus.py`.
- `execution/cmdata/` — the CM library (`sources`, `genes`, `drugs`, `corpus`,
  `qa`, `paths`); `views.py` is retired (views are consumer-owned, ADR 0033).
- `execution/tahoe/build_pseudobulk.py` — Tahoe pseudobulk builder moved from
  `04` (ticket #14): reads `pert_data/Tahoe-100M/source/`, writes the corpus.
- `execution/_parts/` — 14 per-plate CSR parts + `_state/shard_plate.json`;
  local-only resumable cache, not an artifact.
- `execution/build_record/` — tracked build QA/manifests: `tahoe_pseudobulk_*`
  (ticket #14 move) plus the corpus build records tickets #15/#16 add.
- `execution/tests/` — unit plus built-artifact contract tests.
- `output/corpus/v2/` — the one output class: five h5ad corpora, seven tables,
  one coverage README. `tahoe.h5ad` (67,018 × 62,710, 25.8 GB, local-only) is
  here now; tickets #15/#16 build the other four corpora and the tables.

Run from the repository root; code finds the repo root by walking up to
`.git`, never `parents[N]`/`../..` (ADR 0030 §4). Do not build a corpus or
rewrite a ticket #15/#16 product from this skeleton state.
