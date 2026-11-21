# Development comparison of chemical-only Tahoe transfer

[Issue #265](https://github.com/Moke110/arc/issues/265) compares transferring
the TranSiGen ECFP4 encoder alone, transferring the complete Tahoe model,
and initializing all weights randomly.

## Parameter ownership

The native backbone's `feat_embeddings` maps 2,048 ECFP4 bits to 400
features with a linear layer, BatchNorm, ReLU and dropout. Its seven state
tensors contain linear weight and bias, BatchNorm weight and bias, running
mean and variance, and the batch counter. The chemical-only handoff loads
all seven tensors strictly, including three BatchNorm buffers. It leaves
the recipient's control expression encoder, both expression decoders, the
separate dose embedding, and both predictive heads freshly initialized
under the common seed. The predictive heads jointly consume 100 control
latent coordinates, 400 chemical coordinates and 32 dose coordinates.
The chemical encoding therefore has a separate load boundary and was
learned jointly with its source response map.

All model parameters train after initialization in each arm. Each arm
gets a fresh Adam optimizer and a target dose standardizer fitted on
training rows. Only full transfer uses the verified source-to-target
affine dose conversion. Chemical-only transfer has no dose-coordinate
dependency. Exact assessment scores were checked after serializing the
PCC-selected state and loading it into a fresh model for every arm.

## Protocol

Use the same cold-drug-dose fold 0 and three disjoint parent roles as the
[#263 selection audit](optimization-selection.md). Training uses 1,674
profiles from 275 parents. Selection uses 98 profiles from 16 parents.
Assessment uses 94 profiles from 15 parents. Both roles assess doses 1
and 10 µM. The outer test and held strips are excluded from scoring and
training. Evaluate the frozen best parent-equal validation Delta PCC
checkpoint under a common budget of 540 updates (20 epochs) at learning
rate `1e-4`, separately for seeds 42 and 17. Fix native sampling seed
22617 for prediction. Average scores within parent and weight parents
equally. The two seeds reuse the same assessment population.

The nearest-training-parent similarity is maximum Tanimoto over the
frozen 2,048-bit ECFP4 fingerprints of the 275 fit parents. The fixed
bins are low [0, 0.4), medium [0.4, 0.6) and high [0.6, 1]. A bin with
insufficient support is reported with its denominator.

## Independent assessment

All 94 profiles and 15 parents have defined scores in every arm. The random
and full assessments exactly reproduce the matched #263 settings. All three
arms use the same update budget. The selected epochs differ because each
arm selects independently by validation PCC.

| Seed | Initialization | Selected epoch | Δ PCC | Δ RMSE |
| ---: | --- | ---: | ---: | ---: |
| 42 | random | 19 | 0.21897 | 0.46384 |
| 42 | full | 20 | 0.21711 | 0.46642 |
| 42 | chemical-only | 19 | 0.25456 | 0.45250 |
| 17 | random | 20 | 0.21234 | 0.46379 |
| 17 | full | 20 | 0.22668 | 0.46592 |
| 17 | chemical-only | 13 | 0.26105 | 0.45239 |

Chemical-only transfer improved PCC over random initialization by +0.03559
and +0.04872 for seeds 42 and 17, with RMSE changes of −0.01135 and −0.01140.
Compared with full transfer, PCC improved by +0.03745 and +0.03437, and
RMSE fell by 0.01392 and 0.01353. Lower RMSE indicates smaller expression
magnitude errors under this metric.

| Seed | Nearest fit Tanimoto | Parents / profiles | Random Δ PCC / RMSE | Full Δ PCC / RMSE | Chemical-only Δ PCC / RMSE |
| ---: | --- | ---: | ---: | ---: | ---: |
| 42 | low | 9 / 56 | 0.2164 / 0.4618 | 0.2363 / 0.4624 | 0.2575 / 0.4490 |
| 42 | medium | 4 / 26 | 0.2555 / 0.4898 | 0.2172 / 0.4911 | 0.2689 / 0.4876 |
| 42 | high | 2 / 12 | 0.1575 / 0.4211 | 0.1307 / 0.4354 | 0.2126 / 0.3983 |
| 17 | low | 9 / 56 | 0.2092 / 0.4597 | 0.2497 / 0.4614 | 0.2687 / 0.4476 |
| 17 | medium | 4 / 26 | 0.2519 / 0.4833 | 0.2124 / 0.4995 | 0.2665 / 0.4866 |
| 17 | high | 2 / 12 | 0.1472 / 0.4430 | 0.1515 / 0.4192 | 0.2159 / 0.4054 |

Nearest similarities ranged from 0.120 to 0.731. Chemical-only PCC exceeded
both comparators in every bin at both seeds. Its RMSE exceeded the random
arm for the seed-17 medium bin. The high bin represents two parents and
provides limited evidence for highly similar compounds. The pilot supports
retaining this handoff for a longer confirmation across more parent
splits. It leaves open whether the advantage persists beyond 540 updates,
at other doses, with a different source checkpoint, or when copying only the
linear layer without BatchNorm state.

## Reproduction

Use the accepted Step 03 pair for `PAIR`, the pinned native checkout for
`SOURCE`, and the accepted Step 04 experiment for `EXPERIMENT`. Write
to fresh `EVIDENCE` and `REPORT` paths from the repository root in the arc
environment.

```bash
STEP=stage2/04_tsg_transfer
"$STEP/execution/capped_phase.sh" "$EVIDENCE" chemical-development \
  python "$STEP/execution/chemical_transfer.py" \
  --inputs "$PAIR/tenv" --tahoe-inputs "$PAIR/tahoe_bulk" --source "$SOURCE" \
  --experiment "$EXPERIMENT" --output "$REPORT"
```

Run date 2026-11-21, capped CUDA wall time 60.84 s, peak RSS 2,204,276 KiB.
Report SHA-256 `6e16ea6e999cc340920979b5676e46ecbca1a960fe44f4cb82f24d07d92b6ee4`.
The report records training curves, role and manifest digests, input
identities, exact reload checks and similarity denominators.
