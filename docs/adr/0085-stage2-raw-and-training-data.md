# 0085. Publish source-native raw data and plate-matched training data

**Status** — Accepted, 2026-11-24. [Implementation and acceptance](../../stage2/03_data_train/README.md).

Stage 2 training needs one expression row per plate-level perturbation condition. Tahoe and tenv both have conditions that recur on different plates. Keeping a row for each plate and cell line gives every perturbation effect a single experimentally matched reference. This supersedes ADR 0083's active expression publication boundary while preserving its source-data and gene-axis boundaries.

- `01_data_raw` replaces `01_corpus`. Preserve the seven current sources, native gene axes, source observation identities, metadata and integer counts. Retain controls, blanks and unmatched observations. Tahoe bulk retains the existing cell-to-well-and-cell-line aggregation.
- `03_data_train` replaces `03_view`. Publish `output/tenv/` and `output/tahoe-bulk/` on the frozen ordered G_5000 axis. Pool eligible vehicle counts by plate and cell line, using tenv DMSO and Tahoe DMSO_TF. Pool treated counts by plate, cell line, canonical drug parent and measured dose. Project before summation and apply CP10K followed by natural `log1p` to each pooled profile on its selected-axis total.
- Each dataset publishes one control profile vector `C`, a control matrix `M_C` keyed by plate and cell line, and a treatment matrix `M_T` keyed by plate, cell line, parent and dose. `C` is the equal-weight mean of the normalized rows of `M_C`. Each row of `M_T` carries normalized treatment expression, treatment minus matched control, and six fold labels. Keep a membership and exclusion ledger. Exclude incomplete treatment identities, missing compatible controls and groups with zero selected-axis totals from the applicable training artifact.
- Preserve each source's three frozen Step 07 cold-drug and three frozen cold-drug-dose memberships. Map the existing assignments onto the pooled groups only where member assignments agree. Record `train`, `valid` and `test` in the fold columns, and `excluded` in cold-drug-dose where applicable. Keep the held-parent and held-dose exclusion reasons. Study-specific splits remain study-owned.
- Keep existing full-view inputs and sealed records available for replay. Bind the new outputs to the raw generation, gene-axis order, and fold provenance. Verify the candidate before switching the stable path.

See [ADR 0083](0083-stage2-corpus-axis-view-boundaries.md), [ADR 0077](0077-five-k-only-axis-local-normalization.md) and [ADR 0078](0078-dual-source-axis-local-dose-round1.md).

[Implementation specification](https://github.com/Moke110/arc/issues/305).
