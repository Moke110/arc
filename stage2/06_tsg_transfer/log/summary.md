# Accepted Tahoe-to-tenv comparison

Tahoe validation selected fold 2 at epoch 551 (treated-expression MSE 0.02611569).

| tenv fold | Held dose (µM) | Test profiles | tenv-only | Tahoe-only | Transfer | Transfer − tenv-only | Transfer − Tahoe-only |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.2 | 469 | 0.251319 | 0.067867 | 0.287614 | +0.036295 | +0.219746 |
| 1 | 1.0 | 468 | 0.265040 | 0.056841 | 0.271065 | +0.006025 | +0.214223 |
| 2 | 10.0 | 465 | 0.251747 | 0.055017 | 0.224385 | -0.027362 | +0.169368 |
| Equal-fold mean | — | 1402 | 0.256035 | 0.059909 | 0.261021 | +0.004986 | +0.201113 |

Scores are sample-equal Delta PCC within each fold. The final row weights folds equally. Its profile count is the total across the three folds. The [plot table](../figure/data/transfer_scores.csv) records undefined PCC counts by path and fold.

Transfer − tenv-only measures the effect of Tahoe initialization under the fixed tenv training protocol. Transfer − Tahoe-only measures the effect of tenv adaptation from the locked Tahoe checkpoint.

The [acceptance receipt](acceptance.json) binds all nine fits, nine test evaluations and fresh-process replays.
