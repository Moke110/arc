# Compare Tahoe training strategies across three tenv split regimes

**Status:** proposed (2026-11-08). The three split regimes, training seeds
and primary comparison are agreed. The warm-start split unit and training
protocol remain open.

`stage2/05_pdiff_transfer/` asks whether Tahoe pretraining improves
PerturbDiff predictions of tenv drug responses. All test responses come from
tenv. Evaluate `warm-start`, `cold-cm-drug` and `cold-drug`. Within each fold
of a regime, all arms share the frozen tenv split. Each trained model uses
seeds `17`, `42` and `2026` in every fold. Evaluate each deterministic
reference once per split. `no_change` predicts zero Δ and `delta_mean`
predicts the mean of the tenv training group Δ values.

The arms are `no_change`, `delta_mean`, `Tahoe` (Tahoe-only training), `CM`
(tenv-only training), `Transfer` (Tahoe pretraining followed
by tenv fine-tuning with separate stage recipes), and `Mix` (pooled Tahoe and
tenv training in a single fit). Deliver validation-selected checkpoints for the
trained models, learning curves, held-out metrics and PR curves, following the
presentation scope of Step 04.

The two cold regimes use the frozen Step 04 tenv fold 0 train, validation and
test membership and its three `CM` checkpoints. Steps 03 and 04 already
scored these test groups. `cold-cm-drug` holds out drug parent-connectivity
blocks in tenv and permits overlap in Tahoe. `cold-drug` additionally excludes
the tenv validation and test blocks from Tahoe training and validation for all
Tahoe-dependent arms. A warm-start split changes the tenv training membership
and requires its own `CM` fits.

Within each split regime, the primary contrast pairs `Transfer` and `CM`
by fold and seed. Compute the difference in their equally weighted positive
and negative AP@100 and report the individual paired differences and mean.
Retain all eight Step 04 metrics and both directions of the held-out PR curves.
Seed variation describes repeated fits on a fixed split.
