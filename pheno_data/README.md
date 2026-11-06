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

The current consumer is `stage3/01_corpus/` under [ADR 0044](../docs/adr/0044-stage3-source-first-corpora.md).

Heliste is archived with `python pheno_data/Heliste-2021/script/download.py` and
verified offline by `bash pheno_data/Heliste-2021/script/verify.sh`. Dataset S1
contains 689 unique labels and 2,756 supplied relative CTG values per oxygen
condition; Dataset S3 has 48 selected hit labels. Dataset S2 is shRNA data and
is excluded from the chemical archive after bounded acquisition timeouts.
The source assessment is frozen in its declared `local/assessment_v1.json`.
The source-first builder retains resolved chemicals without significance or endpoint
admission. Unknown normalization denominators and independent units remain source
limitations. The previous empty-after-QC assessment is historical evidence only.
HR timing is 24 h pretreatment, 48 h hypoxia and 24 h reoxygenation, with drug
persistence unknown; parallel normoxia is duration matched, about 96 h elapsed.
Frozen chemical resolution is declared in `local/chemical_resolution_v1.json`;
acquisition is separate from the offline `stage3/01_corpus/execution/corpus.py` build.

| Active item | Source evidence | Corpus limitation |
| --- | --- | --- |
| `Sirenko-2013` | Full 24 h beat-count XLS and methods | Parent/use identity, matched control moments and numerical tests unresolved |
| `Heliste-2021` | Full Dataset S1, selected-hit S3, methods and supplement PDF | CTG denominator and independent statistics unresolved; parallel 96 h normoxia isolated |
| `Linn-2024` | Original MTT Figure 3 and article | Exact dose grid and matched vehicle unresolved; 48 approximate summaries retained |
| `Liu-2021` | XML-confirmed Figure 2 and methods | Matched administered vehicle unconfirmed; all eight comparisons excluded |
| `Duan-2013` | Original Figure 2 and article HTML | SYA parent identity and required statistics unresolved; 4 h toxicity panel excluded |

Use `python pheno_data/<Item>/script/download.py` for acquisition and
`bash pheno_data/<Item>/script/verify.sh` for offline verification. Liu and Duan
record use of original primary bytes acquired earlier in the run, rather
than claiming a later refetch succeeded. Figure-derived values retain their
extraction precision and independent review; summary n never expands into wells.

## Huang 2022 — public evidence archive, full screen pending

`Huang-2022` archives the public combined article and supplement, not the full numerical screen. Acquire with `python pheno_data/Huang-2022/script/download.py` (requires `curl` and Python `pymupdf`); verify with `bash pheno_data/Huang-2022/script/verify.sh`. Active registry status means the archived evidence is hash-verified, not that the full dataset or corpus is available.

The eligible baseline imaging design is 24 h at 100 µM in 13 human iPSC-CM donor lines, with Hoechst/PI whole-well live/dead quantification. Public Table S3 supplies selected-hit flags, not continuous full-library values. Full screen acquisition needs an author request; the source's 2,374/2,375 compound discrepancy and reported 10% DMSO remain unresolved. See the [access audit](../docs/research/stage3-huang-2022-access-audit.md), source annotation and unsent request in the item's declared `local/` additions. No baseline release was built.
