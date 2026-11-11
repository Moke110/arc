# Step 05 · Tahoe transfer to tenv

> **Dose encoding limitation.** Tahoe dose embeddings use log-dose interpolation with endpoint clamping. The 10 µM tenv dose receives the learned 5 µM Tahoe row. Transfer starts from this mapping and can update the dose embeddings on tenv. Improving the cross-dataset encoder is a priority for follow-up.

All arms predict held-out tenv drug–dose group means. The cold regimes reuse the Step 04 tenv fold 0 test set, which was also scored in Step 03. Warm-start uses its own held-out parent–dose combinations. Each regime is a fixed fold, so seed variation describes repeated fits on that fold.

## Retrieval reporting

The previous directional-rank comparisons were removed and remain in Git.
Current tenv retrieval uses the three full-CM significance metrics defined by
[ADR 0059](../../../../docs/adr/0059-stage2-significance-retrieval.md).
See the [current release](../../output/retrieval/manifest.json).

## Tahoe parent overlap

Counts reflect the realized Tahoe source memberships and the tenv roles. Mix pools the Tahoe training membership. Tahoe and Transfer use the source training and validation memberships.

| Regime | Tahoe train ∩ tenv validation | Tahoe train ∩ tenv test | Tahoe validation ∩ tenv validation | Tahoe validation ∩ tenv test |
| --- | ---: | ---: | ---: | ---: |
| warm-start | 5 | 9 | 0 | 3 |
| cold-cm-drug | 1 | 2 | 0 | 0 |
| cold-drug | 0 | 0 | 0 | 0 |

## Checkpoint selection

Learning curves and the [checkpoint index](checkpoints.json) retain the original
validation-loss selections. Retrieval rescoring does not select new checkpoints.
