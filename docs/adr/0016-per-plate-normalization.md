# Per-plate matched-control and QN grouping in the rebuilt finetune data

**Status: accepted**

v1's frozen DToxS finetune data (built outside the repo) and the v1.1/v2
rebuilds preprocessed each finetune dataset as **one dataset-wide group**: a
single quantile alignment to the L1000 pretrain distribution and a matched
control computed from *all* control wells of the dataset pooled together. But
every CM finetune dataset is plate-structured with **per-plate control wells**
(dtoxs_primary CTRL/CTRL-DMSO/CTRL-EtOH on plates 1/2/4/5; dtoxs_ipsc CTRL on
plate 6; gse244740 DMSO+NoTreatment per optimization plate; gse262419 VEH+Media
per plate), and the concentration-response assays (gse244740, gse262419) are
explicitly designed so treated and control wells share a plate. Pooling
controls across plates mixes plate batch effects into the control baseline —
the very confound per-plate controls exist to remove.

We decided to make **per-plate normalization** the method in the *rebuilt*
finetune data: each sample's matched control is the mean of the control wells
on its own `(dataset, plate)`, and the 978-landmark matrix is quantile-aligned
to the L1000 pretrain **per `(dataset, plate)` group** (one monotonic transform
per plate). This applies to v1.1, v2, and the new v3 (ADR 0017). **v1's frozen
`input/dtoxs.h5ad` is untouched and remains the dataset-wide historical
anchor** — the difference is recorded as an auditable delta in the step-03
parity/delta reports, exactly like the dose and cell-name deltas before it
(ADR 0015).

Why not keep the DToxS arms dataset-wide to stay v2-identical: every treated
plate carries its own control wells (verified, 0 affected samples), so the
per-plate rule is uniformly applicable and the alternative would be keeping a
known batch-effect confound in the primary/iPSC arms for the sake of
unchanged numbers. The DToxS arms therefore also switch to per-plate; the
rebuilt h5ads' sample *identity* is unchanged (same wells, same order), so
fold columns — generated X-independently from `pert_idx` + seed (ADR 0014) and
imported byte-identically across workflows — survive the rebuild unchanged.

Consequences:

- v1.1's `dtoxs_lincs.h5ad` and v2's `cm_union.h5ad` X/X_ctl values change;
  their archived benchmark numbers become stale and steps 05–11 are re-run.
- Step 02's consensus baseline is **not** affected: it is a single
  dataset-wide reference row for scoring, not a per-sample training control,
  and stays byte-identical to v1's baseline.
- v1 stays frozen in every respect; the agentic report step (`12_report`) is
  deleted from v1.1 and v2 as part of this change (spec 1).
