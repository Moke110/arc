# AGENTS.md — pheno_data

Phenotype external items follow ADR 0030/0038; read `README.md` and `registry.json`
before source acquisition or curation. Active items must verify offline.

- `Sirenko-2013/` — 24 h beat-count workbook and article methods; `author_log.json`.
- `Heliste-2021/` — HL-1 CTG drug Dataset S1, methods and supplementary directory; `author_log.json`.
- `Linn-2024/` — 48 h MTT curves in hiPSC-CM and AC16; `author_log.json`.
- `Liu-2021/` — icariin CCK8 normoxia and H/R; `author_log.json`.
- `Duan-2013/` — primary neonatal-rat MTT after A/R; `author_log.json`.
- Acquisition state in `registry.json` is distinct from corpus admission.
  Originals stay in `source/`; declare each local addition in its author log.

Consumer: `stage3/01_data_corpus/`. Build inputs are frozen; network belongs only
in item acquisition scripts. See ADRs 0038–0041 for corpus scope.
