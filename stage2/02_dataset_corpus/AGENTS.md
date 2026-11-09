# AGENTS.md — stage2/02_dataset_corpus

Routing index. Human record: `README.md`; coverage/provenance/gates:
`output/corpus/v2/README.md`; decisions: ADR 0033, ADR 0026. Status: #14 moved
everything here (`tahoe.h5ad` included), #15 built the four CM corpora, #16 wrote
the seven tables + coverage README, #17 `--verify` **G1–G6 PASS** (G4 reproduces
5/5 frozen v1 content hashes), #18 moved the CM-schema view to `05`; #20 docs.

- `execution/build_corpus.py` — single entry point: `--corpus
  all|dtoxs_primary|dtoxs_ipsc|sctox|tenv`, `--check`, `--repeat-check` (G5),
  `--tables` (seven tables + coverage README + `corpus_coverage.json`),
  `--all`, `--verify` (G1–G6, reads only from disk; record in
  `build_record/verify.json`). `.tmp`→rename, `track_times=False`.
- `--shared-cm build|verify` — ADR 0052 derived corpus at
  `output/corpus/shared-cm-v1/`; bounded native raw-count projection and replay
  in `execution/cmdata/shared_cm.py`, CLI fixtures in `execution/tests/test_shared_cm.py`.
- `execution/cmdata/` — CM library: `paths`, `sources`, `genes`, `drugs`,
  `annotation`, `corpus`, `qa` (G1 + spot check), `corpus_h5ad` (native-axis
  build/digest/write), `tables` (crosswalk/coverage/README, G2/G3). Views are
  consumer-owned; the v1 view writer is frozen in `stage2/_frozen/cm-data-v1/`.
- `execution/tahoe/build_pseudobulk.py` — Tahoe pseudobulk builder (#14);
  `execution/_parts/` — per-plate CSR cache (local-only).
- `execution/build_record/` — tracked records (#14 `tahoe_pseudobulk_*`, #15
  corpus build/reproducibility/exclusions, #16 `corpus_coverage.json`);
  `*.csv` and `output/` are gitignored — force-add.
- `execution/tests/` — `test_cmdata.py`, `test_corpus_contract.py` (G1),
  `test_tables_contract.py` (G2/G3); skip when the local h5ads are absent.
- `output/corpus/v2/` — five h5ad corpora (four CM dense + `tahoe.h5ad` 25.8 GB
  CSR, local-only), `tables/` seven CSVs, `README.md`.

Run from the repository root; find it by walking up to `.git`, never
`parents[N]`/`../..` (ADR 0030 §4). Check `free -g` ≥ 12 GiB and record peak
RSS: a corpus build peaks ~7.7 GiB; `--tables` pays the same global load.
