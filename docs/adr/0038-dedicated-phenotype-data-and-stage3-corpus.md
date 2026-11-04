# Dedicated phenotype provenance category and Stage 3 corpus workspace

**Status: accepted workspace boundary (2026-11-04); corpus partition and sample-grain decisions under discussion**

After completion of the Stage 3 wayfinder map, the user explicitly chose a separate
`pheno_data/` category for cell drug-perturbation phenotype datasets and
`stage3/01_data_corpus/` for quality control, cleaning and standardized corpus
production. This separates experimental phenotype sources from the transcriptomic
perturbation sources consumed by Stage 2.

This decision supersedes the source placement under `pert_data/` and the planned
`01_phenotype_database` step name in ADR 0037. It extends the external-source
category list in ADRs 0001/0030; the per-item provenance layout is retained:
one upstream product per item, `script/` for downloading and offline verification,
untouched `source/`, individually declared `local/` additions where needed,
tracked `author_log.json`, and a source-only `manifest.sha256`.
No existing transcriptomic source is relocated by this boundary decision.

The user requested one corpus file per `(condition, phenotype)` and records of
`(drug, dose, endpoint)`, defining phenotype as the measured indicator and endpoint
as direction (+1/0/−1). Those terms are recorded in the glossary. The precise
condition boundary, whether endpoint belongs to sample identity or is its response,
and admission of unresolved labels into training exports remain live design
questions. No unagreed sample-key or labeling rule is fixed by this record.

This is the grill-with-docs design step. Implementation follows the repository's
spec/ticket workflow after the remaining design decisions are resolved; this record
does not claim production data have been downloaded, labels assigned, or the
qualified 500-observation gate passed.

See [ADR 0001](0001-category-taxonomy.md),
[ADR 0030](0030-external-data-provenance-and-step-layout.md),
[ADR 0037](0037-stage3-phenotype-release-layout.md),
[GLOSSARY.md](../../GLOSSARY.md), and the
[completed wayfinder map](https://github.com/Moke110/arc/issues/21).

## Subsequent corpus-boundary decisions (2026-11-04)

The user resolved fine condition partitions, a common MTT/CCK8 `cell_viability`
phenotype, endpoint-inclusive sample identity, corpus-stage label derivation and
source `readout` retention in
[ADR 0039](0039-condition-corpus-and-model-view-separation.md).
Labeling criteria, same-key repeated measurements and unresolved endpoints remain
under discussion; the workspace boundary above is unchanged.
