# Tahoe representation freeze

A fresh process verified 352 distinct canonical parent structures and all
65,318 treated profiles against the sealed Tahoe input and three cold-drug
folds. All seven arms have a finite native float32 vector for every parent.
The verification found no missing, duplicate, invalid or structurally
mismatched parents. The ECFP4C bits match the existing Tahoe ECFP4 contract.

| Representation | Parents | Width |
|---|---:|---:|
| ECFP4C | 352 | 2048 |
| ECFP6C | 352 | 2048 |
| FCFP4C | 352 | 2048 |
| FCFP6C | 352 | 2048 |
| E3FP | 352 | 1024 |
| KPGT | 352 | 2304 |
| UniMol | 352 | 512 |

E3FP produced one conformer for 14 parents, two for 13 and three for 325.
The capped CPU extraction took 28 min 29.22 s and peaked at 1,993,064 KiB RSS.
The [freeze manifest](representation-freeze.json) binds the vector, parent
and conformer bytes and records the model sources. The [independent
verification](representation-freeze-verification.json) records the complete
fold counts and the exact KPGT matrix digest used for the Step 04 comparison.

Read the [execution guide](../../AGENTS.md#representation-input-freeze) to
reproduce the checkout-local `output/representations/tahoe/` freeze.
