# Heliste 2021 normalization and direction-label audit

Reviewed 2026-11-06 against the archived original article XML, all five pages of
the supplementary-information PDF, Dataset S1 and Dataset S3. This is an evidence
audit, not a labeling-policy change or a rebuilt corpus.

## Primary sources

- [Article, Methods: Drug library screening; Results: Drug library screen; Figure 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC8371130/).
  Archived as `pheno_data/Heliste-2021/source/PMC8371130.1.xml`.
- Supplementary information: `pheno_data/Heliste-2021/source/41598_2021_96033_MOESM1_ESM.pdf`.
  Page 5 describes the datasets; pages 2–4 contain shRNA, zebrafish and Western-blot figures.
- Dataset S1: `pheno_data/Heliste-2021/source/41598_2021_96033_MOESM2_ESM.xlsx`.
  Sheet1, B3:E3 and F3:I3 identify the oxygen arms; B4:I4 identify their four doses;
  rows 5–693 contain the chemical measurements.
- Dataset S3: `pheno_data/Heliste-2021/source/41598_2021_96033_MOESM4_ESM.xlsx`.
  Sheet1 contains selected survival-promoting drug names, not dose-level statistics.

## What normalization is established?

The drug-screen methods explicitly identify 0.01% DMSO-treated cells as vehicle
controls and state that the CTG data were normalized using DMSO-control values.
Thus the reference treatment is known: DMSO vehicle. Calling the screen reference
entirely unknown would be too strong.

The exact transformation is not stated: the inspected sources give no equation,
background correction, control aggregation, plate mapping or explicit statement
of whether each oxygen arm uses its own vehicle denominator or both use a common
normoxic denominator. S1 has values on a roughly percentage-like scale but no
normalization legend defining 100. It has one visible sheet, no cell comments,
no formulas and no vehicle rows. Its title's description as raw viability data
does not make these raw RLU: the article describes normalization and the column
headings explicitly say relative viability.

Figure 2b reports a 22% reduction in vehicle-control viability after H-R. Figure
2c identifies a horizontal line as mean H-R vehicle viability. This is evidence
that the vehicle level depends on oxygen condition in that analysis and cautions
against assuming 100 is the H-R reference. A common normoxic-vehicle scale is a
possible interpretation, not a recovered S1 normalization formula. Neither an
approximate 78 inferred from 22% loss nor a graph-digitized mean is a recovered
plate-matched control for each S1 value.

The supplementary PDF supplies no additional chemical-screen formula or control
table. The separate validation experiment in Figure 4a explicitly rescales using
normoxic and H-R medians; that validation transformation must not be transferred
to the library screen without evidence.

## Can chemical versus vehicle labels be derived now?

Not as a complete, evidence-backed per-dose -1/0/1 target. Missing information
includes numerical matched vehicle references, their variation, replicate units,
and per-dose comparison statistics. S1 supplies one number per chemical, dose
and oxygen arm, without establishing whether it is an individual well or summary.
The four doses are not replicates. Drug-treated normoxia is not vehicle-treated
H-R and cannot substitute for that comparator.

The screen methods select a representative best-performing concentration and
call survival improvement using a value more than one control SD above the H-R
vehicle level. The results also require no reduction in normoxic viability for
the cardioprotective category. This is a selected drug-level screen rule, not
a general per-dose three-class test. Dataset S3 cannot label every dose of a hit
+1 or label every non-hit 0 or -1. The validation ANOVA/t-test/FDR procedures do
not supply missing statistical evidence for the full screen.

Threshold-based labels could be defined after the reference scale is confirmed:
for a justified ratio r = treated / matched vehicle, a possible policy is -1
below 1-delta, 0 within [1-delta, 1+delta], and +1 above 1+delta. This is an
illustrative policy, not an approved rule or a statistical significance claim.
Its 0 would mean within a chosen effect-size band, not demonstrated equivalence.
Using an exact sign with 0 only at equality would make rounding/noise define the
classes. Missing comparison evidence must not silently become class 0.

## Information needed to finish

Ask the authors or screen provider (Misvik Biology) for the S1 normalization
equation; reference oxygen arm and reference value; background subtraction and
plate-control handling; treated and vehicle well-level RLU with plate/oxygen/dose
identifiers; and replicate structure. At minimum, numerical matched reference
values and the transformation are needed for effect thresholds. Replicate data
or adequate summary statistics and the test design are additionally needed for
statistically supported labels. No author message was sent during this audit.

Current corpus measurements remain unchanged. The existing annotation's
oxygen-specific-denominator uncertainty remains justified; its DMSO reference
should be understood as established separately from that uncertainty.
