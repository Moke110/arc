# Stage 2 optimization decision for the next tenv experiment

November 21, 2026 · [Map #260](https://github.com/Moke110/arc/issues/260) · [Decision #266](https://github.com/Moke110/arc/issues/266)

**First confirm chemical-only Tahoe initialization under a longer matched tenv
budget.** Keep random initialization as the deployment reference while that
confirmation is pending. Preserve Mean_cell, audit control reuse and assess
magnitude error alongside PCC. The experiments on this map use frozen G_5000
and development partitions. Published outer-test results enter only as
historical descriptions.

## What the measurements establish

- In published Step 03 cold-drug-dose, Tahoe's profile-weighted latent
  embedding Δ PCC exceeds tenv's by 0.2431. The Mean_cell gap is 0.2408.
  After subtracting each domain's own Mean_cell, the gap in paired increments
  is +0.0022 with a 95% parent-bootstrap interval [−0.0204, +0.0235].
  Parent weighting gives +0.0034 [−0.0192, +0.0242]. In the cold-drug
  scenario, the adjusted profile-weighted gap is larger at +0.0339
  [+0.0164, +0.0499]. Tahoe has more fit profiles and more cell contexts
  while tenv has more fit parents and more wells per condition. The domains
  share no measured dose or cell context in these validation comparisons.
  These results quantify the gap relative to a baseline and leave the effect
  of changing training sample size or assay context unidentified.
  [Coverage audit](../../03_tsg_dose/log/261-coverage/README.md)
- Development repeat PCC rises when the same vehicle pool is subtracted from
  both treated profiles. With same-batch treated pairs and equal-size vehicle
  pools, the increase is +0.1275 to +0.1559 in tenv and +0.3337 to +0.4712
  across eligible Tahoe fold-dose cells. Each Tahoe cell in this construction
  contains only one parent. The tenv data have no natural independent-batch
  repeat pairs. Residual PCC remains positive after subtracting the
  train-only Mean_cell vectors from both prediction and observation.
  Cold-drug-dose residual PCC spans 0.1402–0.2474 in tenv and
  0.3969–0.4621 in Tahoe across validation fold-dose cells. This
  transformation cancels the shared control before computing PCC. The
  coverage audit computes paired differences between PCC scores on
  published outer rows. Residuals can include dose and other effects
  shared across parents. Repeat agreement is a measurement diagnostic,
  and the experiments leave the contribution of shared controls to the
  domain gap unidentified.
  [Control and residual audit](../../03_tsg_dose/log/262-control-reliability/README.md)
- Published Step 04 full transfer gains +0.0363, +0.0060 and −0.0274 Δ PCC
  across the three outer folds, for +0.0050 with equal fold weights.
  The held dose changes with the held parents. In the matched 540-update
  development comparison, full transfer at `1e-5` is about −0.044 versus
  random under either selection rule. Both arms and both seeds choose
  `1e-4` in the equal-budget search. The mean paired PCC difference is then
  +0.0055 under MSE selection and +0.0062 under PCC selection, with the sign
  varying by seed and slightly higher transfer RMSE. Switching selector
  changes the paired PCC difference by at most 0.0015. The historical full
  transfer fits reached lower validation MSE after roughly two to three
  and a half times as many updates as tenv-only. Most pilot choices lie
  at the 20-epoch cap, and only the MSE-best historical weights were
  preserved. The long-budget effect of changing the selector remains open.
  [Accepted Step 04](summary.md) · [Selection audit](optimization-selection.md)
- With the same 15 assessment parents each at 0.2, 1 and 10 µM, the
  full-transfer PCC differences are +0.0183, −0.0038 and −0.0065.
  Every parent-bootstrap interval crosses zero, including the differences
  between endpoint extrapolation and 1 µM interpolation. Full transfer
  increases the fraction of genes with sign errors by about 0.034–0.037
  across doses, and its RMSE is higher at 10 µM by +0.0115
  [0.0005, 0.0227]. Those intervals condition on the two fitted seeds.
  Rotating the held dose also changes which two concentrations are available
  to train the dose mapping. The study does not identify a dose effect
  specific to the source initialization or explain the published fold sign
  changes.
  [Fixed-parent dose crossing](264-dose-crossing/README.md)
- With the same 540 updates and assessment split as the selection audit,
  chemical-only transfer beats random by +0.0356 and +0.0487 Δ PCC
  at seeds 42 and 17 and lowers RMSE by about 0.0114 in both. It also
  beats full transfer by +0.0375 and +0.0344 PCC. The checkpoint copies
  the ECFP4 encoder and its BatchNorm state. Both seeds cover the same
  15 parents and two doses (1 and 10 µM). The nearest-training-parent
  similarity bins contain 9, 4 and 2 parents. PCC gains occur in all three
  bins for both seeds, while chemical-only RMSE worsens against random
  in the medium bin for seed 17. The pilot estimates the effect of that
  exact initialization bundle and budget. It leaves the separate roles of
  learned weights, BatchNorm, and downstream response coupling unresolved.
  [Chemical transfer audit](optimization-chemical-transfer.md)

## Interpretation from the literature

Ben-David et al., Theorem 2, link source-domain classification error,
hypothesis-class divergence and the best joint classifier's error to
target-domain error. This gives a reason to check both feature support and
whether a response rule serves both domains. Tahoe and tenv change measured
dose grids, cell contexts and assay together. Here the classification result
serves as a qualitative guide. We did not estimate its terms for continuous
expression or apply it to Δ PCC.
[Ben-David et al., 2010](https://doi.org/10.1007/s10994-009-5152-4)

Ahlmann-Eltze et al. found that simple mean and linear predictions
were competitive with deep models on their genetic perturbation tasks.
They also observed gains when importing perturbation representations learned
on relevant perturbation data. Their work motivates paired comparisons
against train-only references and testing representation transfer separately.
Their methods distinguish Pearson delta's sensitivity to direction from
an error measure that also captures scale. These results guide our design
for chemical perturbations and provide no tenv effect estimate.
[Ahlmann-Eltze et al., 2025, Fig. 2 and Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12328236/)

## Ranked decision cards

| Priority | Investment | Evidence and trade-off |
| ---: | --- | --- |
| 1 | Confirm the chemical encoder handoff with a common longer budget | The most consistent pilot gain is in chemical-only transfer and includes lower aggregate RMSE. Three matched initializations reuse an existing Tahoe checkpoint and avoid new source fits. Existing assessment has 15 parents and two seeds. |
| 2 | Measure sensitivity to repeat aggregation and control choice | Reuse of vehicle pools materially changes repeat PCC. Run a CPU-only analysis using parent × dose aggregation and equal-size disjoint vehicle pools on matched development conditions. Treat independent-batch tenv reliability as a separate data requirement if the observed controls cannot support it. |
| 3 | Retain both checkpoint selection metrics in the next fits | Recording both metrics has low incremental cost and lets one fit serve both rules. The 540-update pilot favored `1e-4` under both rules. Any further learning-rate grid should use the same candidates and update limits for all initializations and an independent assessment split. |
| 4 | Defer wider dose-architecture changes until the handoff confirms | The fixed-parent crossing did not resolve a dose-dependent full-transfer PCC gain. Some direction and magnitude results favor random initialization. A new crossed three-dose comparison of chemical-only transfer would cost additional fits and should follow the first confirmation. |

## Minimum confirmation and stopping rules

1. Freeze a new parent-disjoint development split inside Step 03 fold 0
   cold-drug-dose development (1 and 10 µM). Keep the original validation
   parents for selection and reserve 30 new assessment parents from the
   frozen training role by a published hash ordering. Require both doses
   for each selected assessment parent. All 275 eligible frozen training
   parents currently cover both doses. Report coverage by dose and Tanimoto
   band. Freeze the actual identities before training. The
   formal test parents, held dose, and excluded strips stay outside this run.
2. Run random, full transfer and chemical-only at `1e-4` with fresh Adam,
   the same 64-profile batch protocol, seeds 42, 17 and 2026, and matched
   checkpoints at budget caps of 540 and 2,160 optimizer updates. Evaluate
   validation at common update steps and retain the earliest maximum
   parent-equal Δ PCC state within each budget. Log the minimum-MSE state
   as a secondary comparison from the same trajectories. Evaluate both
   only after locking their choices. The 2,160-update comparison is
   primary, with nine fits under a single capped GPU lease scheduled
   serially. The historical validation optima show why a gain here
   cannot settle the effect under the original 1,000-epoch budget.
3. Prespecify the primary contrast as chemical-only minus random
   parent-equal Δ PCC at 2,160 updates, with doses weighted equally
   within parent and the three seeds averaged equally. Pair identical
   assessment profiles, report defined denominators and parent-cluster
   95% intervals across the fixed seeds, and show each seed separately.
   Report paired RMSE, full-transfer contrasts, train-only Mean_cell
   differences and residual PCC as secondary outcomes. Recompute the
   comparison with equal well counts per condition and disjoint
   same-batch control pools where the audit shows support. Publish
   missing counts and treat any control sensitivity as descriptive.
4. Advance chemical-only to a formal protocol only if its mean PCC gain
   versus random is at least +0.01, the paired parent-bootstrap lower
   bound exceeds zero, the PCC gain is positive for at least two of three
   seeds, and its mean paired RMSE difference is at most zero. Also require
   its mean PCC to be at least that of the full-transfer arm. These are
   proposed investment thresholds, set before the confirmation. Keep the
   result inconclusive on failed identity checks, fewer than 30 complete
   assessment parents, or a gain that changes sign under the matched
   control sensitivity. If the statistical or magnitude gates fail, stop
   this confirmation after the 2,160-update cap, retain random as the
   reference, and report the finite-budget result. If chemical-only
   fails to establish a gain and most selected states hit the cap, use
   that as evidence that the longer-budget question is still open.
