# Formal corpus statistical families, Parquet export and view-level coverage

**Status: accepted (2026-11-04)**

The user accepted the final three design branches following ADRs 0038–0040.
The formal corpus uses source-aware statistical comparisons, native numerical
lists in Parquet, and coverage reporting without discarding small condition files.

## Statistical route and multiplicity

Prefer a source-reported numerical p-value whose comparison, method, sidedness
and correction status are known and support the intended contrast. Do not apply
BH a second time to an already adjusted value. An author-adjusted value is usable
as the chosen BH result only when its procedure and family match the contract;
otherwise obtain raw p-values or valid inputs for recomputation.

When recomputation is needed and its assumptions are applicable, use a two-sided
Welch test for independent groups or a paired t-test for a documented paired
design. Technical measurements do not create independent experimental units.
Summarization needed for a test does not replace the retained corpus readouts.
Independent-group mean/sample-SD/n can support Welch recomputation; marginal
paired summaries without pairing or difference uncertainty cannot recover a
paired test. Unknown units, insufficient inputs or invalid test results remain
untestable rather than acquiring fabricated p-values.

The BH family is `(source, full condition, phenotype, assay_method,
original experiment batch)`. Include the family's eligible drug/dose tests before
filtering by magnitude or raw-p significance. Record family membership, size,
method and comparison provenance. Selectively reported significant p-values alone
do not establish a complete family. Stars/ns without sufficient numerical inputs
do not become exact p-values. Unavailable significance follows ADR 0040's admission
rule, preserving fully numeric corpus endpoints.

## Export and coverage

Write one Parquet file per `(condition, phenotype)`, with JSON manifest and schema.
Retain source-specific endpoint-bearing sample keys and numerical readout lists,
including per-measurement source/assay/unit/comparison evidence. This replaces
ADR 0037's CSV primary-corpus carrier for the formal design; external original
bytes and their provenance still follow ADR 0030 under `pheno_data`.

Do not impose a 500-row minimum on each fine-condition corpus file. Retain small
files and report coverage. The third-category 500-observation target remains
visible for sources and subsequent model-view selection/aggregation. Report full
sample-key rows alongside drug–dose combinations with source and endpoint removed,
as well as distinct drugs; source or endpoint splitting is not proof of biological
independence. No view is built or qualified at the corpus-design stage.

The design frontier is now resolved. Missing source evidence is handled by the
declared admission/exclusion rules rather than a new unresolved design branch.
The consolidated [formal corpus design](../research/stage3-formal-corpus-design.md)
is the implementation-spec handoff; no production ingestion or spec publication
is performed by accepting this record.

See [ADR 0038](0038-dedicated-phenotype-data-and-stage3-corpus.md),
[ADR 0039](0039-condition-corpus-and-model-view-separation.md),
[ADR 0040](0040-numeric-endpoints-and-source-readout-sequences.md), and
[GLOSSARY.md](../../GLOSSARY.md).
Statistical API definitions: [Welch](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html),
[paired test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_rel.html),
[summary statistics](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind_from_stats.html),
[BH](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html).
List storage: [Apache Arrow Parquet type handling](https://arrow.apache.org/docs/python/parquet/parquet_type_handling.html).
