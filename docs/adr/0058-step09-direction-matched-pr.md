# Match significance direction and stop PR at attainable recall

**Status:** accepted (2026-11-11). Amends the significance-retrieval definition
in ADR 0057 following clarification of the intended metric.

A Step 09 true positive must be significant at the frozen DESeq2 adjusted
p-value threshold and match the sign of its DESeq2 log2 fold change. A wrong
sign is a false positive and leaves the true response unrecalled. Rank
nonzero predictions by absolute group Δ; zero predictions abstain. Tied
magnitudes enter together. Retain all significant genes in the recall
denominator, including genes outside coverage and genes predicted with the
wrong sign. `@full` uses the CM axis; `@G` uses the model axis.

For displayed PR, pool selected and correctly matched gene–group pairs at a
common absolute-Δ threshold across a fixed set of groups with defined AP.
Precision is pooled hits divided by pooled predictions, and recall is pooled
hits divided by the fixed total significant-gene count. Stop at the last
attainable threshold; do not extrapolate to recall 1 or remove groups as
recall rises. This avoids the changing-population artifact of the earlier
reach-conditioned macro display. Even a full-axis model can stop below recall
1 because wrong directions are unrecoverable. Raw precision need not be
monotone. Keep exact per-group AP and its equal-group macro mean in the
metric tables; the displayed pooled PR area is a different statistic.

The Round 1 test set was already examined. Preserve the original validation
promotion and unsigned score/figure release as historical evidence. Rescore
the same frozen predictions for validation and test under the revised metric,
without rerunning model selection or treating revised validation comparisons
as new pre-test decisions. Revised test labels still bind the historical seal.

## See

- [Original ablation contract](0057-step09-xpert-landmark-reference.md)
- [Step 09 release and reproduction](../../stage2/09_xpert_optimize/README.md)
