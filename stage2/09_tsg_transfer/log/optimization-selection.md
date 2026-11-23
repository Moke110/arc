# Development audit of validation selection and fine-tuning budget

[Issue #263](https://github.com/Moke110/arc/issues/263) compares selection by
treated-expression MSE with selection by parent-equal Delta PCC and repeats
the comparison under a small, equal-budget learning-rate search.

## Preserved Step 04 evidence

The accepted fit and validation-replay receipts verify against the frozen
tenv cold-drug-dose split. All six target fits used seed 42, Adam at
`1e-5`, and 27 optimizer updates per epoch. The original validation role
uses held-out parents within the two training doses in each fold.

| Fold | Arm | Validation profiles | Validation parents | Epoch 1 MSE | Best MSE | Best epoch | Updates at best | Total epochs |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | tenv-only | 192 | 31 | 0.587980 | 0.215537 | 231 | 6237 | 281 |
| 0 | transfer | 192 | 31 | 0.874985 | 0.204846 | 566 | 15282 | 616 |
| 1 | tenv-only | 186 | 31 | 0.589427 | 0.224279 | 230 | 6210 | 280 |
| 1 | transfer | 186 | 31 | 0.982381 | 0.209906 | 809 | 21843 | 859 |
| 2 | tenv-only | 188 | 31 | 0.555878 | 0.176057 | 310 | 8370 | 360 |
| 2 | transfer | 188 | 31 | 0.824367 | 0.167144 | 632 | 17064 | 682 |

Transfer started with higher validation MSE and reached lower validation MSE
in all three folds. It took 2.0 to 3.5 times as many updates to reach its
own minimum. The histories contain only native training loss, treated MSE
and the improvement/step record. Each fit preserves one validation-best
checkpoint. The Delta PCC-optimal epoch and same-trajectory reselection
effect for these accepted runs cannot be recovered from those artifacts.

## Bounded development comparison

The experiment fixes fold 0 development rows and excludes the outer test
and both outer held strips. It uses 1 and 10 µM. Training uses the
original 1,674 profiles from 275 parents. The 31 frozen validation parents
are ordered by SHA-256 of parent InChIKey and alternated between selection
(98 profiles, 16 parents) and independent assessment (94 profiles, 15
parents). Each role retains its frozen profile order. The split manifest
SHA-256 is
`6a716fdb129dcad124c626a398cf0afc58f57b85ba35453ec51d2ba28cf1ccf0`.
The source-only Tahoe winner checkpoint SHA-256 is
`76c9f5846e06932309ece2801d204511b0b278314b6f3afb25d99708d5fd1b0e`.

Both initializations use fresh native Adam, the same two learning rates
(`1e-5`, `1e-4`), seeds 42 and 17, and exactly 540 optimizer updates per
fit (20 epochs). The target dose transform is fitted on training rows and
the transfer arm uses the verified affine dose bridge. Validation prediction
uses sampling seed 22617. Each selection rule uses the earliest epoch in case
of a tie. Treated-expression MSE weights profiles equally. PCC and RMSE
average profiles within parent and then give each defined parent equal
weight. Both criteria use the same validation predictions. Assessment is
performed after restoring each selected state in memory. It uses independent
parents from that selection and never drives epoch or learning-rate choice.

### Fixed-recipe comparisons on independent assessment

Pairs are tenv-only / transfer. Differences are transfer minus tenv-only.
Each selected epoch corresponds to 27 times that many optimizer updates.
All 94 assessment profiles and 15 assessment parents have defined PCC in
every row.

| Seed | Learning rate | Selection | Epoch pair | Delta PCC pair | PCC difference | Delta RMSE pair | RMSE difference |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 42 | 1e-5 | MSE | 20 / 20 | 0.1253 / 0.0879 | -0.0374 | 0.5339 / 0.7544 | +0.2205 |
| 42 | 1e-5 | Delta PCC | 19 / 20 | 0.1249 / 0.0879 | -0.0370 | 0.5399 / 0.7544 | +0.2145 |
| 42 | 1e-4 | MSE | 18 / 20 | 0.2205 / 0.2171 | -0.0034 | 0.4647 / 0.4664 | +0.0017 |
| 42 | 1e-4 | Delta PCC | 19 / 20 | 0.2190 / 0.2171 | -0.0019 | 0.4638 / 0.4664 | +0.0026 |
| 17 | 1e-5 | MSE | 20 / 20 | 0.1382 / 0.0882 | -0.0500 | 0.5344 / 0.7559 | +0.2215 |
| 17 | 1e-5 | Delta PCC | 20 / 18 | 0.1382 / 0.0875 | -0.0507 | 0.5344 / 0.7649 | +0.2305 |
| 17 | 1e-4 | MSE | 20 / 20 | 0.2123 / 0.2267 | +0.0143 | 0.4638 / 0.4659 | +0.0021 |
| 17 | 1e-4 | Delta PCC | 20 / 20 | 0.2123 / 0.2267 | +0.0143 | 0.4638 / 0.4659 | +0.0021 |

At the formal learning rate of `1e-5`, the 540-update transfer arm
underperformed tenv-only by 0.0437 Delta PCC with MSE selection and 0.0438
with Delta PCC selection, averaged over the two paired seeds. The two
rules chose different epochs for three of the eight fits. Across the four
same-rate and same-seed comparisons, reselection changed the paired PCC
difference by at most 0.0015 in absolute terms.

### Equal-budget learning-rate search

For each initialization and seed, compare the two learning rates using
the same criterion used to select the epoch. The independent assessment
set contributes no values to that choice.

| Seed | Selection | Learning-rate pair | PCC difference | RMSE difference |
| ---: | --- | ---: | ---: | ---: |
| 42 | MSE | 1e-4 / 1e-4 | -0.0034 | +0.0017 |
| 42 | Delta PCC | 1e-4 / 1e-4 | -0.0019 | +0.0026 |
| 17 | MSE | 1e-4 / 1e-4 | +0.0143 | +0.0021 |
| 17 | Delta PCC | 1e-4 / 1e-4 | +0.0143 | +0.0021 |

Both rules selected `1e-4` for both initializations and seeds. The
mean paired PCC differences were +0.0055 (MSE) and +0.0062 (Delta PCC),
with opposite directions across seeds. The same search left transfer
RMSE slightly higher for both seeds. Under this short development budget,
learning rate explains more of the observed contrast than which metric
selects an epoch. Most epochs were selected at the cap, and the historical
MSE optima took thousands of updates. These measurements leave the
long-budget epoch-selection effect unresolved. The paired seeds measure
initialization sensitivity on the same 15 assessment parents.

## Reproduction

Run this command from the repository root in the `arc` environment,
using the accepted Step 03 pair for `PAIR`, the pinned native source for
`SOURCE`, and the accepted Step 04 `experiment` directory for
`EXPERIMENT`. `EVIDENCE` and `REPORT` are new result destinations.

```bash
STEP=stage2/04_tsg_transfer
"$STEP/execution/capped_phase.sh" "$EVIDENCE" selection-development \
  python "$STEP/execution/selection_audit.py" \
  --inputs "$PAIR/tenv" --tahoe-inputs "$PAIR/tahoe_bulk" --source "$SOURCE" \
  --experiment "$EXPERIMENT" --output "$REPORT" --epochs 20
```

The CUDA prototype completed on 2026-11-21 in 72.17 seconds with peak RSS
2,134,580 KiB. The replayable JSON report has SHA-256
`45c246fad024e17b75238b3c4525adb1a84052d7b11be0276c2b953e3e5f408b`.
It contains per-epoch selector MSE, PCC and RMSE, fitted epochs and updates,
role population digests, and the independent assessment of every restored
selected state.

The review replay used the shared canonical Delta vector scorer. It kept
all identities, selected epochs and counts identical, and the largest
absolute numeric difference across the JSON was `8.33e-17`. The replay
took 83.13 seconds at 2,135,100 KiB peak RSS and produced SHA-256
`b343b61d471ed98fe031287588b9bea0c221bfff747c7e713ec7d4bc97915101`.
