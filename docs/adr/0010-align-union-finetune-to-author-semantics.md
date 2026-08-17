# Align the union-CM finetune data to the author L1000 semantics (xpert-fda-miri-v2)

The v2 union-CM finetune (ADR 0009) benchmarked worse than v1 (union warm-start
MSE_deg 7.24 vs v1 1.12, Pearson_deg 0.225 vs 0.321). This ADR records the
**data-alignment mitigation**: the union finetune data is rebuilt so every
sample's dose/time index semantics match the XPert author's L1000 conventions,
each dataset's expression distribution is quantile-aligned to the pretrain, the
finetune config's bin range is reset to the pretrain's, and the retrained model
is benchmarked **apples-to-apples on v1's DToxS test folds**.

## Context

The v2 union mixed five `pert_data/cm` datasets with incompatible semantics and
a dominant platform offset:

1. **Time** — the union renumbered time as `{24:0, 48:1, 144:2}`; the author
   pretrain and v1 use `{3:0, 6:1, 24:2}` (24 h → index 2). The 48 h
   (`dtoxs_primary`, `dtoxs_ipsc`) and 144 h (`gse247541`) samples have no
   author time index at all, and the union's 24 h samples carried a time token
   whose pretrain meaning was "3 h".
2. **Dose** — `dtoxs_primary` recorded a 0.0 placeholder dose, so all its
   samples landed in `pert_dose_idx 0` regardless of any real exposure.
3. **Platform drift + bin scale** — per-dataset CPM left the TempO-Seq majority
   (`gse262419`, ~55% of union samples) ~2 log2 units below the DToxS-family
   mean, and the global bin range had been widened to `[0, 18.166]` (pretrain:
   `[0, 15.002]`) solely to cover TempO-Seq outliers, coarsening and
   misaligning the bins from the pretrained bin embeddings.

## Decision

- **Author time map** — the union is restricted to `AUTHOR_TIME_MAP = {3:0,
  6:1, 24:2}` (`_lib/utils.py`; replaces the step-owned `{24:0, 48:1, 144:2}`).
  **v1 DToxS slot-24 convention** (follow-up decision): the DToxS-family
  datasets are real 48 h experiments, but v1's frozen dtoxs.h5ad encoded them
  at the author 24 h slot, so `dtoxs_primary` and `dtoxs_ipsc` are re-added
  under that same convention (`pert_time 24.0`, `time_idx 2`) instead of being
  dropped for their 48 h exposure. `gse247541` (144 h, 12 treated) exits the
  union entirely. Union composition: `gse262419` (3,783) + `gse244740`
  (2,424) + `dtoxs_primary` (400) + `dtoxs_ipsc` (296) = **6,903 samples**, all
  24 h → `pert_time_idx 2`, identical to v1's DToxS encoding. The screen still
  predicts at 24 h; step 07's predict `--time-idx` default is updated to 2
  (the library `predict.py` already used the author map).
- **No placeholder dose, with the v1 single-dose exemption** — treated samples
  with `dose <= 0` and no real concentration are excluded (`gse262419`'s 32
  "Blank" wells). Exception: `dtoxs_primary` is the v1 single-dose experiment
  (no recorded concentration) and enters at `pert_dose 0 → dose_idx 0`,
  exactly as v1's dtoxs.h5ad did. Every other union sample's `pert_dose_idx`
  is derived from its real µM concentration via the author
  `assign_pert_dose_corrected` bins (0–9).
- **Per-group quantile alignment (QN)** — within each matched-control group
  (plate for `gse262419`/`gse244740`, dataset-wide for the DToxS datasets),
  the group's pooled `log2(CPM+1)` distribution is rank-mapped onto the pooled
  empirical quantiles of the L1000 pretrain h5ad, using one monotonic transform
  for every sample of the group so the trt/ctl contrast is preserved. Applied
  after CPM/log2 and **before** matched-control computation, so `obsm["X_ctl"]`
  is the control mean of aligned values. The reference (a row-subsample of the
  pretrain matrix) is cached to `03_cm_finetune_data/output/qn_reference.npy`.
- **Config reset** — the finetune/benchmark config `max_value`/`min_value` are
  reset to the pretrain values `[0, 15.002]`; after QN the union X lies in that
  range, so bin width (~0.119) and bin semantics match the pretrained bin
  embeddings.
- **Auditability** — the step-03 manifest gains `n_union`, `n_excl_time`,
  `n_excl_nodose` (12 time-out + 32 no-dose classified samples) and the
  per-dataset summary records the dropped datasets; the schema gate is extended
  into the **data-semantics gate** (author time/dose validity, zero dose only
  in the `dtoxs_primary` single-dose convention, X within the pretrain range,
  union composition = the four 24 h-slot datasets, manifest `n_union` matches
  union rows) and aborts on failure.
- **Apples-to-apples benchmark** — the retrained union warm-start fold models
  are evaluated with `mode=test` on v1's frozen `input/dtoxs.h5ad` using its own
  `warm-start_fold-1..10` test masks (test-only config
  `05_xpert_finetune_bench/config/cm_eval_dtoxs.yaml`, union `cm` semantics +
  `[0, 15.002]` bins). Because the aligned union's time encoding equals v1's,
  no re-encoding of the eval data is needed. Per-fold metrics are compared to
  v1's archived warm-start benchmark (`benchmark_dtoxs_folds_vs_v1.csv`).

## Consequences

- **In-distribution (union) warm-start benchmark** — MSE_deg 7.24 → 1.07,
  RMSE_deg 2.69 → 1.03, MAE_deg 1.83 → 0.75, Pearson_deg 0.225 → 0.277,
  Spearman_deg 0.219 → 0.260, R2_deg 0.057 → 0.059, with fold std collapsing
  from ~0.10 to ~0.01. The platform-offset distortion is gone; the union model
  now performs at v1's error scale on its own data.
- **Apples-to-apples (v1 DToxS folds)** — on the identical test set the union
  model does **not** match v1: MSE_deg 1.28 vs 1.12, Pearson_deg 0.086 vs
  0.321, Spearman_deg 0.075 vs 0.181, R2_deg −0.22 vs −0.009. Re-adding the
  DToxS family at the 24 h slot narrowed the gap (from Pearson 0.056 without
  DToxS data), but the union model — in which the DToxS family is only ~10% of
  samples — still generalizes worse to DToxS folds than the DToxS-only v1
  model. This is the honest reading the benchmark was designed to surface, and
  it is the open question for whether the union finetune is preferable to v1's.
- The DToxS family is back in the union under v1's slot-24 convention; only the
  144 h dataset (`gse247541`) is excluded from training. The 48 h/144 h source
  times remain the real values in `pert_data/cm` (read-only) and the raw
  extractors; the slot-24 normalization is a step-03 semantic decision.
- Steps 05/06/07 were re-run (GPU-fp-nondeterministic across environments, as
  before); the archived predictions from the pre-alignment model are stale and
  must be regenerated before re-running 08–12.
- `gse262419`'s platform-missing WRB/FAM57A (zeroed) map into the low tail of
  the aligned distribution, as documented in step 03.
