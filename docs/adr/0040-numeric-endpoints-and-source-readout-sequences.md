# Numeric endpoints and source-aware readout sequences in the formal corpus

**Status: accepted rule structure and grouping (2026-11-04); cutoffs and untestable-input policy under discussion**

The user resolved the next grill-with-docs frontier: drug time is pre-HR treatment
duration for HR experiments; endpoint is assigned numerically using a conjunction
of effect-magnitude and significance criteria; source dataset becomes part of
sample identity; matching records within one dataset retain all readouts as a
sequence. Final corpus endpoint values are restricted to +1, 0, and −1.

## Condition time semantics

`condition_drug_time` is the pre-HR treatment duration when HR follows dosing, and
the direct treatment duration when there is no HR. Hypoxia and reoxygenation
durations are separate `condition_*` attributes. A drug schedule attribute records
continued treatment, washout, redosing or unknown status during those phases.
Known absence of HR has zero hypoxia/reoxygenation durations; missing timing
information remains unknown rather than becoming zero. This does not equate
experimental elapsed time with actual drug exposure.

## Corpus identity and measurement retention

One file still represents one `(condition, phenotype)`. Within that file, the
unique sample key is `(source, drug, dose, endpoint)`, where `source` identifies
the upstream dataset, not each file, plate, replicate or publication of it.
The full identity therefore includes condition and phenotype as well.

Within a dataset, records with all these fields equal are merged into one sample
whose `readout` is a numerical sequence retaining all source measurements. Do not
replace the sequence with its mean. Preserve per-value assay method, units,
representation, original experiment/replicate and source locations so that MTT
and CCK8 can share `cell_viability` without erasing how each value was obtained.
Across datasets, source remains different and samples stay separate. Different
endpoints also stay separate. Later model views choose further aggregation.

Endpoint evaluation must precede the final endpoint-keyed merge. Statistical
comparison groups and the measurements supporting each test remain traceable;
grouping on the outcome must not determine which observations enter that test.
Published summary values do not expand into invented replicate measurements.
The primary treatment of raw versus author-normalized source representations is
still to be specified; no extra biological replicate is created by retaining two
representations of the same measurement.

## Numeric endpoint rule

For an evaluable comparison, assign +1 when the upward magnitude requirement and
the selected p-value significance criterion both pass; assign −1 when the downward
magnitude requirement and significance both pass. Otherwise assign 0. The user
requires all admitted corpus records to have numeric endpoints, with no null or
textual unknown endpoint.

This operational class 0 replaces the earlier requirement for positive equivalence
evidence in ADR 0036 and the wayfinder direction policy. It means the joint change
criterion was not met, rather than proving the exact absence of biological change.
Keep effect and test evidence so downstream users can interpret the label.

Magnitude scale/cutoffs, significance cutoff and raw-versus-adjusted p choice,
test selection/provenance, and handling of unevaluable inputs remain live design
questions. A missing p-value is not a tested nonsignificant p-value. If magnitude
is known to fail, the conjunction is already false; if magnitude passes but the
p-value cannot be established, the input policy must be explicitly decided.
No missing p, control distribution or uncertainty is fabricated by this record.
There will be no nonnumeric endpoint in the final corpus whichever admission
policy is selected.

See [ADR 0039](0039-condition-corpus-and-model-view-separation.md),
[ADR 0036](0036-stage3-continuous-evidence-before-direction-labels.md),
[GLOSSARY.md](../../GLOSSARY.md), and the
[source-evidence audit](../research/stage3-sample-count-audit.md).
