# AGENTS.md — pheno_data

Phenotype external items follow ADR 0030/0038; read `README.md` and `registry.json`
before source acquisition or curation. Active items must verify offline.

- `Sirenko-2013/` — 24 h beat-count workbook and article methods; `author_log.json`.
- Other candidates are recorded in `registry.json`; acquisition state is distinct
  from corpus admission. Each acquired item keeps upstream bytes in `source/`
  and declares every generated addition under `local/` in its author log.

Consumer: `stage3/01_data_corpus/`. Build inputs are frozen; network belongs only
in item acquisition scripts. See ADRs 0038–0041 for corpus scope.
