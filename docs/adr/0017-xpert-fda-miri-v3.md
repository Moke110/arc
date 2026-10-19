# xpert-fda-miri-v3 — all-CM finetune sibling (per-plate normalized)

**Status: accepted**

The finalized MIRI screen (v1, ADR 0007) froze a DToxS-only finetune; its
siblings v1.1 (DToxS-LINCS testbed, ADR 0012) and v2 (primary+iPSC union, ADR
0009/0015) progressively widened the finetune corpus. v3 is the **all-CM
finetune sibling**: a full 11-step workflow under
`stage1/xpert-fda-miri-v3/` built from v2 with the finetune-data scope
extended to **four CM datasets** — `dtoxs_primary`, `dtoxs_ipsc`, `gse244740`,
`gse262419` — under the **per-plate normalization** of ADR 0016. `gse247541`
is excluded (6-day exposure outside the author time map, 2 in-HG drugs, 12
treated samples). Everything else is v2-identical: the screened 1,004-drug
scope, step-02 baseline, MIRI scoring/filters/integration, and the predict
exposure (dose_idx 0 / time_idx 2). v2 stays authoritative and untouched.

Two deliberate deviations from v2 worth recording:

1. **Benchmark: the new datasets are train-only.** v3's benchmark test sets
   stay **byte-identical to v2's** (primary rows import v1.1's fold columns
   verbatim; iPSC rows per the ADR 0015 isolation rules). gse244740/gse262419
   rows never appear in a test set; in warm-start folds they are all-train, and
   in cold-drug folds their rows for that fold's test drugs are excluded
   (`''`) across all four datasets. The benchmark therefore answers precisely:
   *holding the test sets fixed, does adding the two dose-response datasets to
   the training corpus change primary/iPSC test performance?* Column-by-column
   comparable with v1.1/v2.
2. **No agentic report.** v3 is an **11-step workflow** (`01_drug_embed` …
   `11_integrate`); the agentic bilingual `12_report` stage is dropped, and
   spec 1 deletes it from v1.1/v2 as well. v1's frozen `11_report` stays.

The finetune corpus is a large platform mix. The built union is **~4,318 rows
/ ~180 unique drug_idx**: dtoxs_primary 662, dtoxs_ipsc 296, gse244740 ~2,424
(38 drugs), gse262419 ~936 (101 in-HG drugs of 464 chemicals — the TempO-Seq
panel is mostly environmental/toxicology compounds absent from the author
table, so fewer in-HG rows than the initial ~3,500 estimate). gse244740, not
gse262419, dominates the union (~56% of training rows); that platform
imbalance is a documented property of the experiment, not a defect — the
research question is whether the broader corpus helps. The union h5ad is
schema-identical to v1's `input/dtoxs.h5ad` in the fields the XPert loader
consumes, with the same `cm` dataset branch (`num_pert_dose 10`) and
`cm_union.h5ad` naming as v2. GPU steps re-run per spec 2; the final screen
still predicts at the v1 exposure.
