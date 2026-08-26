# Align v2 finetune to v1.1 — primary + iPSC-CM, real-dose bins, imported folds

**Status: accepted**

v2's union finetune (ADR 0009) trained on four CM datasets; the finalized v2 screen
benchmarked self-contained folds that did not match v1.1's. We decided to realign v2
to v1.1: the XPert finetune data narrows to `dtoxs_primary` + `dtoxs_ipsc`
(958 samples / 46 unique drugs), each dataset is preprocessed independently
(978 landmarks → CPM → log2 → QN to the L1000 pretrain → matched control → author
dose binning) and concatenated into `cm_union.h5ad` with a single-cell encoding
(`cell_iname="cardiomyocyte"`, `cell_idx=0`), and the benchmark test sets become
byte-identical to v1.1's fold test sets by importing v1.1's fold columns verbatim
for the primary arm (iPSC rows assigned per isolation rules: warm-start all-train,
cold-drug test drugs excluded from both datasets' train).

Two deliberate deviations from v1/v1.1 worth recording:

1. **Primary dose is now real-concentration binned** (not the v1 single-dose
   `dose 0 → dose_idx 0` convention). Concentrations come from the DToxS-primary
   paper Table S1 (23 KIs), the iPSC-CM series matrix (shared non-KIs), and
   literature-curated clinical Cmax for the 8 primary-only drugs, each with a
   citation (SB-431542 uses its standard in-vitro assay concentration). This
   **amends the "no external database" provenance rule for the dose field only** —
   concentrations must be cited. Consequence: under the author µM bins, 17/46
   primary drugs still land in bin 0 (clinical Cmax < 0.21 µM), so the real-dose
   scheme chiefly separates the higher-exposure drugs (bins 1–9).
2. **Single-cell encoding** is retained (num_cell_id 1) even though two cell
   systems now train together — the model never distinguishes donors/lines.
   Fold isolation gates assert well identity (no identical well in train+test)
   and drug identity (cold-drug test drugs excluded from both datasets).

The final screen still predicts at `dose_idx 0 / time_idx 2` (v1 convention),
and step 03 now audits only the two finetune datasets (gse244740/gse262419/
gse247541 are no longer extracted in the build; the standalone drug-overlap
report still audits all five for the data lake).
