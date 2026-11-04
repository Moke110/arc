# Stage 3 — phenotype direction corpus

Stage 3 consumes separately archived experimental phenotype sources from
[`pheno_data`](../pheno_data/README.md). Its first step,
[`01_data_corpus`](01_data_corpus/README.md), builds one typed Parquet partition
per complete condition and phenotype, retaining source-specific measurements,
comparison evidence and numeric endpoint labels under ADR 0038–0041.

Source acquisition, numerical corpus validity, subsequent model-view coverage and
training readiness are separate outcomes. Model views, Stage 2 prediction
matching and model training belong to later work.
