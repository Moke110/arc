# Compare CPST loss terms on shared tenv evaluation folds

**Status** Accepted design, 2026-11-28.

Rename `06_tsg_shape` to `06_tsg_module`. Replace the active alpha step
with `07_cpst_loss` and expose the common CPST model through the shared
utility API. This supersedes the active Step 07 interface ownership in
[ADR 0099](0099-carps-vae-alpha-and-chemical-representations.md). Keep the
[historical alpha release and replay route](https://github.com/Moke110/carps/tree/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/07_cpst_alpha)
at its pinned Git revision. Preserve sealed products and source identities
for replay, including the Step 06.4 checkpoints under the original layout.

## Evaluation protocols

Publish a versioned tenv fold and seed API in `cpst/util`. Steps 04–09 use
five frozen cold-drug folds and training seeds 17, 42, 43, 2026 and 2027.
Use the tenv design from [ADR 0107](0107-cpst-tahoe-module-transfer.md) so
Step 11 can consume the same ledger. Deterministically balance five groups of
canonical drug parents. Keep every dose, plate and eligible condition for each
parent in its assigned group. Rotate each group through test and validation
once, using the other three groups for training. Match parent, condition and
partition identities across the eight tenv axis views. Freeze split membership
before training, with its own seed and digests. Keep existing Step 03 inputs
and frozen historical fold identities available for replay.

Step 10 uses a separate nine-fold cold-dose-drug protocol through the same
API, with training seeds 17, 42 and 2026. Cross the three frozen tenv
held-parent groups of 153 parents each with 0.2, 1 and 10 µM. For each pair,
test only the intersection of held parents and held dose. Keep the original
training and validation parent assignments for that parent group and use
only the other two doses. Exclude the two strips matching exactly one
holdout. Each dose is tested on all three disjoint parent groups, and every
eligible condition enters test once across the nine folds. This extends
[ADR 0103](0103-carps-vae-dose-inputs.md) and separates the effects of the
held dose and the held-parent group. Report 27 fold-and-seed results per arm
and retain the parent-group and dose identities in the report.

## Step 07 loss comparison

Use Step 06.4's `XΔ-XΔ` backbone, axis-local Denoise-DEG 5k, matched
same-plate DMSO, frozen KPGT and log-dose standardized on the training fold.
Use a single zero-weight reference and cross the four metrics below with
α ∈ {0.1, 0.5, 1, 1.75}. All 17 arms cover 25 fold-and-seed coordinates
each, for 425 registered coordinates.

`L = 0.1 K + MSE_X1 + MSE_Δrec + MSE_Δpred + α Lmetric(Δpred, Δ)`

`K` sums the control and measured-Δ prior KLs and the predicted-to-measured-Δ
KL. Use the same measured control to define observed and predicted signed Δ.
The base MSE terms sum over conditions and genes. At α = 0 use only the
base loss. Use the Step 05.1 scale for the additional term.

| Added term | Reduction |
| --- | --- |
| MSE | Sum squared errors over conditions and genes. The effective prediction MSE weight is `1 + α`. |
| L1 | Sum absolute errors over conditions and genes. At α = 1.75 this is the incumbent CPST loss. |
| Cosine | Sum `1 − cosine` across conditions and multiply by the number of genes. |
| Pearson | Center each predicted and observed response across genes, then apply the cosine reduction. |

Use the Step 05.1 similarity normalization with epsilon `1e-8`. A zero
vector after any centering has similarity zero. Select the checkpoint by the
minimum validation value of the arm's full loss. Label training and validation
learning curves `Loss` throughout the rerun. Pair Step 07 results against α = 0
on the same fold and seed.

## Rerun and release

Rerun all formal arms in Steps 04, 05.1–05.5, 06.1–06.4 (including the
06.1 width extension), 08, 09 and 10. Keep each study's fixed model
configuration, arm definitions, comparison reference and downstream model
choices. The Step 07 loss grid adds evidence while the rerun uses the accepted
baseline for each successor step. Record the comparison metric priorities
in the [CPST README](../../cpst/README.md#reporting-for-the-tenv-rerun).

Register the full roster before fitting. Reuse fits only after matching the
scientific and runtime identity, actual partition membership and seed, and
verifying replay. Recalibrate distribution metrics on the applicable
training partition for each new fold where the study uses calibration.
Independently verify held-out scores and replay published figures twice
from saved CSVs. Promote the verified successor releases and retain a pinned
route for historical results.
