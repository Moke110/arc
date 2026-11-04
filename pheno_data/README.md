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

Heliste is archived with `python pheno_data/Heliste-2021/script/download.py` and
verified offline by `bash pheno_data/Heliste-2021/script/verify.sh`. Dataset S1
contains 689 unique labels and 2,756 supplied relative CTG values per oxygen
condition; Dataset S3 has 48 selected hit labels. Dataset S2 is shRNA data and
is excluded from the chemical archive after bounded acquisition timeouts.
The source assessment is frozen in its declared `local/assessment_v1.json`.
Unknown normalization denominators, independent replicate units and canonical
identities produce explicit scientific exclusions. HR timing is 24 h
pretreatment, 48 h hypoxia and 24 h reoxygenation, with drug persistence unknown.
The parallel normoxia experiment has about 96 h elapsed time. The adapter is
`stage3/01_data_corpus/execution/heliste_adapter.py --output BUNDLE.json`; it retains
all 5,512 values and two recognized CTG partitions without inventing labels.
