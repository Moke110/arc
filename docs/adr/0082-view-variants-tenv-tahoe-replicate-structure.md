# Test whether within-condition replicate structure explains the tenv–Tahoe gap

**Status** — accepted (2026-11-22).

## Context

ADR 0078 sealed the dual-source cold-drug-dose comparison on the paired axis-local
G_5000 views. Under the native `no-dose` TranSiGen arm, tenv reaches Δ PCC 0.2655 and
Δ RMSE 0.421, while Tahoe reaches 0.4225 and 0.155 — a Δ PCC gap of 0.157 and a Δ RMSE gap
of 0.27. That gap is the largest unexplained difference between the two Stage 2 sources.

The sources differ in how many profiles share one condition, where a condition is
`(cell_line, drug parent InChIKey, dose)`. Tenv has 4,206 treated wells over 1,377
conditions with no singleton condition at all: 1,320 conditions carry three same-plate
wells, 48 carry four (each plate's `Replicate_Number = 4` reference compound at all three
doses), and 9 carry six (one parent reported as two salt/acid test articles). Its
`matched_control_group` is one same-plate 12-well DMSO group. Tahoe has 65,318 treated
profiles over 52,777 conditions, and 81.5% of those conditions are single-profile; its
repeats mostly sit on different plates and share a two-well `DMSO_TF` control per plate
and cell line. Within-condition Δ repeat correlation has median 0.121 on tenv and 0.207 on
Tahoe, and the condition-internal Δ noise floor has median 0.343 on tenv and 0.070 on
Tahoe. The sealed Δ RMSE of each source sits essentially on its own noise floor.

A heuristic single-sample Δ PCC ceiling — the profile-weighted mean of the repeat
correlation's square root over replicated conditions — is about 0.38 for tenv and 0.43–0.52
for Tahoe. Tahoe therefore performs near its ceiling and tenv has headroom. Two things
follow. Replicate structure can explain the Δ RMSE gap and part, not all, of the Δ PCC gap.
And a hypothesis that repeated tenv wells "confuse" the model, so Tahoe's mostly
single-profile conditions predict better, is testable within each source but cannot on its
own account for a gap tenv could not close even with perfect training labels.

## Decision

Create the study step `stage2/tsg_view`. It owns an expression-view builder that starts
from the Bases integer corpus and the ADR 0077 axis-local recipe, the TranSiGen fits on
each view, and the report. It is a non-numbered study step: it produces comparison
evidence, holds no downstream weight release, and changes no product of Steps 01–06.

The builder publishes eight view products on the frozen G_5000 axis, all under the same
ADR 0078 cold-drug-dose folds (split seed 226, training seed 42, tenv doses 0.2/1/10 µM
and Tahoe doses 0.05/0.5/5 µM). Every view keeps the source's treated profiles and
matched-control pairing except where named.

- `full` is the ADR 0077 recipe unchanged. It is published as a digest-equal alias of the
  accepted `01_bases` release, rebuilt in staging from corpus integers and asserted
  against that release; the step stores one authoritative expression payload, not a copy.
- `merge` carries one replicate-pooled profile per condition. Pooling sums the retained
  G_5000 integer counts of the member profiles and renormalizes the sum on its own
  selected-axis total; normalized values are never averaged. Its control is pooled the
  same way from the union of the member profiles' matched control groups. Replicates on
  different plates are pooled together, so tenv has 1,377 and Tahoe 52,777 profiles.
- `solo` carries one treated profile per condition, drawn once per source by a seed frozen
  before any fit and independent of folds; its control stays the full-view control. A
  sensitivity product `solo-deepest` repeats the draw by taking the largest
  `library_total` member instead.
- `full-pooled-label` carries the full view's profiles and inputs with each treated
  profile's Δ replaced by its condition's pooled Δ, so label noise can move independently
  of profile count.

The model is native TranSiGen under the ADR 0078 `no-dose` arm
(`tsg-dose-no-dose-native-f32-v1`), with the Step 04 training recipe, validation-MSE
checkpoint rule and sample-equal Δ evaluation. The matrix is the seven main views by three
folds by training seeds 42 and 43, plus three `solo-deepest` fits at seed 42 — 45 fits.

Evaluation truth is frozen across views, so a view changes only its training rows. The
primary truth is each treated profile's Δ against its own matched control, as in ADR 0078.
The secondary truth is the condition-mean Δ of a replicated test condition, computed by
the same replicate-pooling recipe from corpus integers; it removes test-side replicate
noise, and its population is small and reported with its own counts (153 tenv test
conditions; 595 replicated Tahoe test conditions in the inspected fold). Δ PCC, Δ Spearman,
Δ RMSE and Δ R² are reported with equal-fold mean and fold spread, plus per-dose panels
where the held dose is the fold's omitted dose.

The study is pre-registered as a partial-contribution test. Support is a tenv Δ PCC gain
from `merge` or `full-pooled-label` of at most about 0.12 that stays below the single-sample
ceiling near 0.38, with the tenv–Tahoe gap still present. A larger gain, or a gain that
crosses the ceiling, revises the ceiling estimate rather than confirming the hypothesis.
The Tahoe arms are a manipulation control, since 81.5% of Tahoe conditions have a single
profile and are unchanged by pooling or drawing.

Out of scope: dose-method comparison, gene-axis change, edits to Steps 01–06 products,
weight publication, and any attempt to attribute the residual gap to cell-line diversity,
corpus size or assay properties, which vary together across the two sources.

## Consequences

`01_bases` stays authoritative for the full view; the step adds a declared merge
generation whose control recipe differs from ADR 0077, and that difference is recorded in
its manifest. Derived views store expression and control only and derive Δ, so the new
payload is roughly five gigabytes, most of it the two Tahoe derived views. The 45 fits cost
about one hour of one capped GPU scope. The report separates what replicate structure
explains from what it does not, and the heuristic ceiling becomes a named quantity that a
later measurement-error model can replace.

## See

- [ADR 0077](0077-five-k-only-axis-local-normalization.md)
- [ADR 0078](0078-dual-source-axis-local-dose-round1.md)
- [Bases](../../stage2/01_bases/README.md)
- [Glossary](../../GLOSSARY.md)
- [Execution constraints](../agents/execution-environment.md)
