# Preserve fine conditions, endpoint-bearing samples and source readouts before model views

**Status: accepted corpus boundaries (2026-11-04); labeling criteria and repeated-measurement storage under discussion**

The user refined the formal Stage 3 corpus design after ADR 0038. The corpus is
partitioned into one file per `(condition, phenotype)`. Condition attributes use
the `condition_*` field naming convention, including `condition_cell_type`,
`condition_drug_time`, `condition_hypoxia_time`, and
`condition_reoxygen_time`, with further experimentally meaningful context retained.
Exact filename encoding and time-field semantics are not fixed by this record.

Conditions remain fine at the corpus stage. A subsequent model data view can
select and combine conditions and decide how to aggregate multiple samples.
The corpus does not preempt those choices by averaging away original observations.
This replaces using the broader wayfinder counting stratum as the corpus file
partition; the original evidence-counting conclusions remain historical facts.

The user explicitly chose a shared phenotype `cell_viability` for MTT and CCK8.
They may therefore occur in the same condition/phenotype corpus partition. Retain
the assay method, units, source and measurement representation so source readouts
are not silently numerically interchanged. This decision does not automatically
add CTG, Calcein imaging, or other assays to that shared phenotype.

Within a condition/phenotype partition, retain `(drug, dose, endpoint)` as the
sample identity, including endpoint (+1/0/−1). Opposite endpoints for the same
drug/dose are retained as distinct corpus samples; any conflict resolution or
sample merging belongs to the model view. Repeated source measurements with an
identical full sample key must remain traceable; their exact grouping/record
representation is still being discussed, rather than silently averaged.

`01_data_corpus` owns automatic derivation of endpoint from continuous source
measurements under rules to be frozen in this design. The corpus also retains the
actual source measurement as `readout`, with its measurement evidence. This refines
ADR 0036's earlier deferral of label generation to a later effort; it does not
approve any numerical threshold, assume absent variance, or equate nonsignificance
with no change. Handling an unresolved endpoint remains a live design decision.

These decisions refine the corresponding early wayfinder choices in ADRs
0035–0037 and the prior inclusion/counting/direction/curation documents. In
particular, their rules excluding endpoint from corpus identity or separating MTT
and CCK8 as different corpus phenotypes are not the formal corpus design.
No claim about a qualified 500-sample dataset is introduced; the appropriate
counting/release gate for the refined corpus and future view remains to be clarified.

See [ADR 0038](0038-dedicated-phenotype-data-and-stage3-corpus.md),
[ADR 0036](0036-stage3-continuous-evidence-before-direction-labels.md),
[GLOSSARY.md](../../GLOSSARY.md), and the
[original wayfinder curation contract](../research/stage3-curation-contract.md).
