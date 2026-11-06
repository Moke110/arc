# Build Stage 3 corpus samples from source measurements without endpoint gates

**Status: accepted scope changes (2026-11-06); identity edge cases under discussion**

The user reopened the Stage 3 corpus design to retain source measurements without requiring significance tests or endpoint labels. This amends the endpoint-based admission and identity rules in ADRs 0039–0041; source measurements must not be excluded merely because those derived results cannot be obtained.

The requested sample key is `(source, rep, condition, chemical, dose, readout)`: `rep` denotes replication number; `condition` is a state such as baseline or hypoxia–reoxygenation, with more detailed `condition_*` attributes outside this key; `chemical` is null for controls and canonical SMILES for treatments; `dose` is null for controls and a micromolar number for treatments; `readout` is the outcome value supplied by the source. Missing source chemical structures require LLM-assisted handling; the resolution and validation policy is not yet decided.

Each source receives an annotation describing experimental design, control type, readout meaning and other interpretation context. This trades the previous label-ready substrate for a measurement-preserving corpus whose later consumers can choose their own contrasts and labels.

The key is recorded as requested, not yet proven collision-free. Replicate scope, multiple assays or times, readout equality, control subtypes, summary-only evidence, chemical ambiguity and physical export partitions remain to be clarified. Existing code and historical design documents have not yet been migrated; this record does not authorize skipping the subsequent spec and ticket workflow.

See [ADR 0039](0039-condition-corpus-and-model-view-separation.md), [ADR 0040](0040-numeric-endpoints-and-source-readout-sequences.md), [ADR 0041](0041-formal-corpus-statistics-export-and-coverage.md), and [GLOSSARY.md](../../GLOSSARY.md).
