# Stage 3 — source-first phenotype corpora

[`01_corpus`](01_corpus/README.md) builds separate Heliste 2021 baseline and
hypoxia–reoxygenation cardiomyocyte `cell_viability` corpora from frozen
[`pheno_data`](../pheno_data/README.md) sources and chemical-resolution evidence.
The source-first contract is [ADR 0044](../docs/adr/0044-stage3-source-first-corpora.md).
The corpus includes derived ±20% magnitude direction labels under an explicit
same-arm vehicle=100 assumption. [`02_base_morgan`](02_base_morgan/README.md) implements the structure-and-dose
baseline under [ADR 0045](../docs/adr/0045-stage3-morgan-direction-baseline.md).
Its current public workflow covers cold-drug preparation and individual saved
runs; the complete three-setting benchmark remains subsequent work.
