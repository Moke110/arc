# Build Stage 3 corpus samples from source measurements without endpoint gates

**Status: accepted decisions through the second interview round (2026-11-06); remaining representation choices under discussion**

The user reopened the Stage 3 corpus design to retain source measurements without requiring significance tests or endpoint labels. This amends the endpoint-based admission and identity rules in ADRs 0039–0041; source measurements must not be excluded merely because those derived results cannot be obtained.

The sample fields are attributes, not a unique key: `source`, `rep`, `condition`, `chemical`, `dose` and `readout`. No uniqueness constraint is inferred from any combination of them. `condition` is a state such as baseline or hypoxia–reoxygenation; detailed `condition_*` attributes remain available. Controls have null chemical and dose; treatments use canonical SMILES and a micromolar dose.

`rep` is an index distinguishing real, physically separate samples whose other attributes, except readout, are exactly the same. It does not assert biological independence. Physically separate samples remain separate even when their readouts are equal; attribute equality must not trigger merging. A second representation or repeated extraction of the same measurement does not establish another physical sample.

`readout` is a dictionary from readout type to supplied value, for example `individual`, `grouped_mean`, `grouped_sd`, `grouped_n` and `grouped_p`, with further types possible. Preserve supplied grouped evidence without expanding a reported n into invented individual samples. How group summaries are linked to physical samples remains to be decided.

Missing chemical structures use LLM-assisted identity resolution with recorded evidence and deterministic structure canonicalization. Ambiguous identities, mixtures and unresolved stereochemistry enter a review queue; null chemical is reserved for controls. Original source names and resolution evidence are retained.

Each source receives an annotation with an experimental-design overview and row-linked mappings for assay, units, scale/normalization, control subtype, replicate hierarchy, timing and source location. This accommodates heterogeneous experiments within a source and makes readouts interpretable without guessing.

These decisions trade the previous label-ready substrate for a measurement-preserving corpus whose later consumers choose contrasts and labels. Group-summary placement, multiple readouts, physical export partitions and additional representation details remain under discussion. Existing code and historical design documents have not yet been migrated; this record does not authorize skipping the subsequent spec and ticket workflow.

See [ADR 0039](0039-condition-corpus-and-model-view-separation.md), [ADR 0040](0040-numeric-endpoints-and-source-readout-sequences.md), [ADR 0041](0041-formal-corpus-statistics-export-and-coverage.md), and [GLOSSARY.md](../../GLOSSARY.md).
