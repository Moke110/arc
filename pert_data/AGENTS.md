# AGENTS.md — pert_data

Index of the CM drug-perturbation data lake. Eight external 项 (ADR 0030): each holds
`script/` (downloader + offline verifier), `source/` (untouched download, gitignored),
`local/` (declared local additions, gitignored), and `author_log.json` (tracked provenance).

- `Tahoe-100M/` — HF `tahoebio/Tahoe-100M` full atlas, rev `349e2a55…`, CC0-1.0; `author_log.json`
- `LINCS-LDS-1593/` — DToxS primary-CM LINCS LDS-1593 well-level package (manual clue.io download); `author_log.json`
- `DToxS-iPSC/` — DToxS iPSC-CM GEO SuperSeries GSE217424 (GSE217421 Conv + GSE217423 DGE); `author_log.json`
- `GSE244740/` — Rosell-Hidalgo ScreenSeq hiPSC-CM cardiotoxicity panel; `author_log.json`
- `GSE247541/` — Xu antiviral hiPSC-CM panel; ADR 0072 author-positive TOX subset; `author_log.json`
- `GSE262419/` — Tsai TempO-Seq hiPSC-CM environmental-chemical panel; `author_log.json`
- `XPert-L1000/` — XPert L1000 landmark-gene pretrain h5ads (upstream provenance unknown); `author_log.json`
- `Annotation/` — pinned HGNC gene set + public-registry compound snapshot + frozen CM drug audit seed; `author_log.json`

Layout contract: ADR 0030. Category docs: `README.md`. Paths in tracked code point at
`<Item>/source/…` (downloaded bytes) or `<Item>/local/…` (local additions); item names are
final. Every item carries `manifest.sha256` (453.0 GB hashed, 27/27 verified);
`script/verify.sh` checks it offline.
