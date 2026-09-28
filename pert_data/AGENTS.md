# AGENTS.md — pert_data

Index of the CM drug-perturbation data lake (raw/processed transcriptomic datasets + LINCS L1000 pretrain corpus). One line per entry; details live in each dataset's own files/reports.

- `dtoxs_primary/` — DToxS primary CM, released LINCS LDS-1593 well-level package (`Data/`, `Metadata/`, `ProcessingPipeline/`, `ReadMe.txt`, `LDS-1593.tar.gz`, `DToxS_summary.md`) plus frozen `v1-frozen/`.
- `dtoxs_ipsc/` — DToxS iPSC-CM (Hansen 2024): GSE217421 Conv + GSE217423 DGE GEO files (SuperSeries GSE217424).
- `gse244740/` — Rosell-Hidalgo ScreenSeq: processed counts + series matrix.
- `gse247541/` — Xu COVID-19 antivirals: raw counts + series matrix.
- `gse262419/` — Tsai TempO-Seq: 16 per-plate CSVs + hash CSV + series matrix.
- `xpert-l1000/` — XPert L1000 pretrain corpus (`data/` h5ads + `L1000_summary.md`).
- `cm_drugs.html` — pert_data CM datasets in-HG drug overlap report.
- Removed (provenance only): GSE221595 (8 TKIs), GSE146095/GSE146096 raw GEO deposit.
