# Stage 3 uses existing provenance items and inspectable phenotype releases

**Status: accepted (2026-11-04)**

Stage 3 phenotype sources are drug-perturbation experiments. They fit the existing
`pert_data` category and one-upstream-product-per-item provenance contract, without
adding a fifth category or changing what existing Stage 2 consumers read. Original
upstream bytes stay in the item's source tree; declared local extractions and
cross-source releases have separate ownership and checksums.

The future first step is `stage3/01_phenotype_database/`. It adopts ADR 0030's
AGENTS/README/execution/output structure, with human coverage, attrition and gate
tables directly at step level. Immutable releases contain typed UTF-8 CSV tables
and JSON schema, rules and release records. CSV uses RFC 4180 field quoting with
explicit LF record separators, an escaped null sentinel, and preserved source
tokens. This favors inspectable small experimental tables over a database service
or a transcriptome-specific container.

Measurements, representations, experimental combinations, comparisons and evidence
links remain separate. Validation checks provenance, schema, source fidelity,
identity/conditions, deduplication/counting, the third-category 500 threshold, and
contrast/label evidence. Passing a validator means its rules were enforced; it
does not make every retained row qualified. Inventory and honest pilots can exist
while qualified counts or label readiness remain unverified or fail.

The first release does not require author contact. Missing evidence is isolated;
insufficient eligible data does not authorize scope expansion or a lower threshold.
Those are decided failure exits, not unresolved wayfinding branches. Production
ingestion, publication of a spec, and final model/label calibration are subsequent
work, not actions performed by accepting this record.

Detailed fields, nullable/control semantics, derivation inputs, gate definitions and
workflow handoff are in the [curation contract](../research/stage3-curation-contract.md)
and its [owning ticket](https://github.com/Moke110/arc/issues/27).
See [ADR 0001](0001-category-taxonomy.md),
[ADR 0030](0030-external-data-provenance-and-step-layout.md),
[ADR 0035](0035-stage3-stratified-phenotype-curation.md),
[ADR 0036](0036-stage3-continuous-evidence-before-direction-labels.md), and
[wayfinder map](https://github.com/Moke110/arc/issues/21).

## Subsequent workspace amendment (2026-11-04)

The user's subsequent grill-with-docs instruction replaces this record's placement
under `pert_data/` and planned first-step name with a dedicated `pheno_data/`
category and `stage3/01_data_corpus/`. See
[ADR 0038](0038-dedicated-phenotype-data-and-stage3-corpus.md).
The original record remains the historical wayfinder decision; corpus partition,
sample identity and training-label admission are being refined in the new design.

The completed formal design in ADRs
[0039](0039-condition-corpus-and-model-view-separation.md)–[0041](0041-formal-corpus-statistics-export-and-coverage.md)
also replaces the CSV primary carrier, endpoint-free sample identity, nullable
corpus labels and per-stratum 500 gate from the early wayfinder plan.
See the current [formal corpus design](../research/stage3-formal-corpus-design.md).
