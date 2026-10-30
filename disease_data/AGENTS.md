# AGENTS.md — disease_data

Shared per-GEO MIRI disease inputs. Three external 项 (ADR 0030): each holds `script/`
(the R recipe + offline verifier), `source/` (upstream GEO-derived files, gitignored),
`local/deg/` (locally derived DEG lists), and `author_log.json` (tracked provenance).

- `GSE108940/` — rat MIRI (Agilent GE2-v5); limma WT IR vs Sham; rat→human ortholog mapping; `author_log.json`
- `GSE137293/` — pig infarct scar vs border zone (no sham control); Welch t-test; pig→human mapping; `author_log.json`
- `GSE160516/` — mouse MIRI (Affymetrix Clariom S); limma pooled IR vs Sham; primary rank reference; `author_log.json`

Each 项 carries its own copy of the shared `script/{map_to_human.R,map_orthologs.R}` helpers
(no cross-项 dependency) plus its series-specific `script/build_deg.R` and `script/verify.sh`.
`local/deg/` keeps the tracked per-GSE DEG JSONs and the regenerable `*_l1000.csv` tables
(declared in `local_additions[]`; the global `*.csv` rule keeps those out of git).
Every item carries `manifest.sha256` (442 GB hashed, 27/27 verified); `script/verify.sh`
checks it offline. Category overview: `README.md`. Layout contract: ADR 0030.
