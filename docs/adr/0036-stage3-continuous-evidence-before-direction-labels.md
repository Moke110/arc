# Stage 3 retains continuous evidence before assigning direction labels

**Status: accepted (2026-11-04)**

Public phenotype tables often omit vehicle values, experimental uncertainty, or
normalization semantics. Treating those gaps or nonsignificant results as the
no-change class would embed an unsupported target into future model training.

The database therefore retains measurements, author-normalized representations,
continuous contrasts when justified, and their evidence. The direction label may
remain unresolved. Unknown is not class 0 or a fourth model class; the intended
model output remains increase/no change/decrease. Readout direction also does not
define therapeutic benefit.

Every contrast must identify a matched biological/injury and vehicle context.
A number near 100 does not establish a reference or formula. Missing raw values,
assay zeros, no-peak outputs, and source-derived zeros remain distinguishable.
MTT, CCK8, CTG, beat count, and other endpoints are not silently pooled.

Final numerical direction/equivalence margins and interval algorithms belong to a
subsequent controlled labeling/model specification. Until its frozen rule and
evidence requirements are satisfied, the default direction is unresolved. This
permits an honest data foundation without manufacturing balanced or complete
three-class labels from incomplete experiments. Qualified sample counts and
label-ready counts are separate gates.

The detailed policy and cases live in the
[direction-evidence decision](../research/stage3-direction-evidence-decision.md)
and [owning decision ticket](https://github.com/Moke110/arc/issues/26).
See [ADR 0034](0034-stage3-phenotype-data-boundary.md),
[ADR 0035](0035-stage3-stratified-phenotype-curation.md), and
[GLOSSARY.md](../../GLOSSARY.md).

## Subsequent corpus-design refinement (2026-11-04)

The user subsequently assigned endpoint derivation to `01_data_corpus` and unified
MTT/CCK8 under the corpus phenotype `cell_viability`; see
[ADR 0039](0039-condition-corpus-and-model-view-separation.md).
The labeling criteria and unresolved-endpoint policy are being refined in the
formal design. No threshold or no-change inference is approved merely by moving
label generation into the corpus step.

The user's next decision in [ADR 0040](0040-numeric-endpoints-and-source-readout-sequences.md)
requires numeric corpus endpoints and assigns 0 when an evaluable comparison
fails the joint magnitude-and-significance rule. This replaces this record's
equivalence-only class-0 and nullable-corpus policies for the formal corpus; handling
of missing test evidence is being specified explicitly rather than inferred.
