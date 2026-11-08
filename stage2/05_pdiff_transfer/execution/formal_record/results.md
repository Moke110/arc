# Step 05 · Tahoe transfer to tenv

> **Dose encoding limitation.** Tahoe dose embeddings use log-dose interpolation with endpoint clamping. The 10 µM tenv dose receives the learned 5 µM Tahoe row. Transfer starts from this mapping and can update the dose embeddings on tenv. Improving the cross-dataset encoder is a priority for follow-up.

All arms predict held-out tenv drug–dose group means. The cold regimes reuse the Step 04 tenv fold 0 test set, which was also scored in Step 03. Warm-start uses its own held-out parent–dose combinations. Each regime is a fixed fold, so seed variation describes repeated fits on that fold.

## Transfer minus CM

The paired difference uses the equal mean of positive and negative AP@100. A regime mean requires all three seed pairs. Group-level scores and undefined reasons are in the full scoring output. The compact release records [undefined reason counts](undefined_reasons.json) and each metric’s valid-group counts in [per_run.csv](per_run.csv). The cold CM and deterministic references share the `cold-shared` scope.

| Regime | Learned slots | Seed 17 | Seed 42 | Seed 2026 | Complete pairs | Mean difference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| warm-start | 12/12 | -0.01737 | -0.03511 | -0.03219 | 3/3 | -0.02822 |
| cold-cm-drug | 12/12 | -0.02825 | -0.02096 | -0.01864 | 3/3 | -0.02262 |
| cold-drug | 12/12 | -0.02966 | -0.02252 | -0.02481 | 3/3 | -0.02566 |

## Tahoe parent overlap

Counts reflect the realized Tahoe source memberships and the tenv roles. Mix pools the Tahoe training membership. Tahoe and Transfer use the source training and validation memberships.

| Regime | Tahoe train ∩ tenv validation | Tahoe train ∩ tenv test | Tahoe validation ∩ tenv validation | Tahoe validation ∩ tenv test |
| --- | ---: | ---: | ---: | ---: |
| warm-start | 5 | 9 | 0 | 3 |
| cold-cm-drug | 1 | 2 | 0 | 0 |
| cold-drug | 0 | 0 | 0 | 0 |

## Results and checkpoint selection

- [Eight-metric summary](summary.csv) and [held-out results](metrics.png)
- [Directional precision–recall plots](pr_curves.png) and [per-run sources](pr_points.json)
- [Learning curves](learning_curves.png) and [checkpoint index](checkpoints.json)
- [Publication checksums](publication.json)
