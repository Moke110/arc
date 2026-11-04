# Phenotype source archive

`registry.json` is the phenotype item inventory. Active means acquired and subject
to mandatory offline integrity checks. Candidate and unavailable records preserve
acquisition gaps; they are never counted as verified empty datasets.

Each item follows [ADR 0030](../docs/adr/0030-external-data-provenance-and-step-layout.md):
`script/` acquisition and verification, untouched ignored `source/`, individually
declared `local/` additions, tracked `author_log.json` and source-only
`manifest.sha256`. A source archive does not imply corpus eligibility.

Acquire Sirenko with `python pheno_data/Sirenko-2013/script/download.py` and verify
with `bash pheno_data/Sirenko-2013/script/verify.sh`. The public workbook has
131 compound entries and four 24 h doses; its Control column is a library
classification, not vehicle measurements. Article rights are reserved; these
locally downloaded bytes are not redistributed in git.

The consumer is `stage3/01_data_corpus/`; the frozen statistical contract is in
[the formal corpus design](../docs/research/stage3-formal-corpus-design.md).
