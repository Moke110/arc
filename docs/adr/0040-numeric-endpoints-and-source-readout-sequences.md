# Numeric endpoints and source-aware readout sequences in the formal corpus

**Status: accepted rule structure, grouping and parameterization (2026-11-04); statistical-family and export details under discussion**

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

## Accepted parameterization (2026-11-04)

The user accepted all four proposed execution choices in the next interview round:

- Fold change is relative to the matched-control mean. The initial operational
  magnitude boundaries are `FC >= 1.2` upward and `FC <= 0.8` downward. These are
  project starting parameters, not claimed to be calibrated biological boundaries.
- Significance uses Benjamini–Hochberg-adjusted `p_adj < 0.05`. Preserve raw p-values
  and adjusted values with their method, comparison and family provenance. Do not
  substitute an unadjusted threshold without an explicit rule change.
- Exclude a candidate from the corpus when its magnitude cannot be established,
  or when magnitude passes but the required p-value criterion cannot be evaluated.
  Retain source material and the exclusion reason. Missing p is not silently 0,
  1, or a nonsignificant result. Admitted corpus endpoints remain numeric.
- Prefer raw source values for the `readout` sequence. Preserve supplied normalized
  values separately. If only normalized values exist, retain those source values
  as readout and mark their representation; do not invent raw values. Distinct
  representations of one measurement do not create extra repeats.

For an evaluable comparison, this gives the rule:

```text
FC >= 1.2 and p_adj < 0.05  -> endpoint = +1
FC <= 0.8 and p_adj < 0.05  -> endpoint = -1
otherwise                  -> endpoint = 0
```

The last branch applies only when the conjunction is decidable. If FC is known
to lie strictly between 0.8 and 1.2, magnitude already fails and endpoint is 0
without inventing a p-value; record unavailable test evidence as unavailable.
If FC crosses a magnitude boundary but the significance criterion is unavailable,
exclude the candidate rather than executing a programming-language NaN comparison
and falling through to 0. Nonpositive/invalid control denominators do not yield
an invented fold change.

The p-value computation/acceptance route, BH family scope, handling of incomplete
published test families, final file format and the placement of the historical
500-observation coverage gate remain the next design frontier. No particular test
or export format is implicitly approved by accepting the numbers above.

The chosen BH cutoff is a labeling procedure, not a claim that the final
amplitude-filtered label subset has exactly 5% false discovery rate. The numerical
implementation must preserve procedure and assumptions; official definitions are
in [SciPy's false-discovery-control documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html).
