# AGENTS.md — stage3/01_data_corpus

- `README.md` — public workflow and corpus contract.
- `execution/` — adapters, offline build/verify, tests and execution evidence.
- `output/` — immutable versioned condition/phenotype Parquet releases.
- Source changes belong in `../../pheno_data/`; preserve original bytes and
  register each local curation addition under ADR 0030.
- Statistics or sample-key changes: read ADR 0039–0041 and
  `../../docs/research/stage3-formal-corpus-design.md` first.
