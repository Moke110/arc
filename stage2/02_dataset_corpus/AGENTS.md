# AGENTS.md — stage2/02_dataset_corpus

Routing index. Human record: `README.md`; contracts: ADR 0033 (the corpus
layer) and ADR 0026 (CM identity). Status: ticket #15 built the four CM
corpora; #16 tables, #17 `--verify`, #18 downstream view, #20 docs are open.

- `execution/build_corpus.py` — single entry point: `--corpus
  all|dtoxs_primary|dtoxs_ipsc|sctox|tenv`, `--check` (raw source counts),
  `--repeat-check` (gate G5), `--verify` (TODO, ticket #17). One corpus at a
  time, `.tmp`→rename, `track_times=False`.
- `execution/cmdata/` — CM library: `paths` (layout + corpus contract counts),
  `sources` (four raw loaders; `cell_line`/`n_raw_rows` provenance), `genes`,
  `drugs`, `annotation`, `corpus`, `qa` (G1 gates + source spot check),
  `corpus_h5ad` (native-axis build/digest/write). Views are consumer-owned: the
  retired v1 view writer is archived in `stage2/_frozen/cm-data-v1/execution/`.
- `execution/tahoe/build_pseudobulk.py` — Tahoe pseudobulk builder (moved from
  `04` in ticket #14).
- `execution/_parts/` — 14 per-plate CSR parts + `_state/shard_plate.json`;
  local-only resumable cache, not an artifact.
- `execution/build_record/` — tracked records: `tahoe_pseudobulk_*` (#14),
  `corpus_build.json` + `excluded_compounds.csv` + `corpus_reproducibility.json`
  (#15). CSV is gitignored, force-add it.
- `execution/tests/` — `test_cmdata.py` (units) and
  `test_corpus_contract.py` (G1 + counts; skips when local h5ads are absent).
- `output/corpus/v2/` — the one output class: five h5ad corpora (four CM dense
  + `tahoe.h5ad` 25.8 GB CSR, local-only), seven tables (#16), one README.

Run from the repository root; code finds the repo root by walking up to `.git`,
never `parents[N]`/`../..` (ADR 0030 §4). Check `free -g` ≥ 12 GiB, build one
corpus at a time (peak measured ~7.7 GiB), record peak RSS.
