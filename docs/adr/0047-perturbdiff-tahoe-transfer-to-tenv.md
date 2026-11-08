# Compare Tahoe training strategies across three tenv split regimes

> **Important limitation to address in follow-up: cross-dataset dose encoding.**
> Tahoe's fitted dose levels (0.05, 0.5 and 5 µM) and tenv's (0.2, 1 and
> 10 µM) have no overlap. Step 05 temporarily bridges the learned categorical
> dose embeddings by linear interpolation in log(µM), clamped to the nearest
> endpoint outside the source range. Thus the zero-shot `Tahoe` arm receives
> the 5 µM embedding for a 10 µM tenv case. `Transfer` starts from the same
> mapping and can update it with tenv data. Interpret the zero-shot scores as
> conditional on this mapping, and record the limitation in the Step 05 README
> and results. Improving dose encoding across datasets is a priority for
> follow-up work.

**Status:** proposed (2026-11-08). The three split regimes, training seeds
and primary comparison are agreed. Warm-start validation membership and the
Mix training recipe still need to be frozen.

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
Tahoe-dependent arms. Use one separate `warm-start` fold, with one held-out
dose per tenv drug parent and test groups balanced across dose levels. Keep
all records for a parent–dose combination together, including aliases and
repeated wells. Every validation and test parent remains represented in
training at another dose. Fit three new `CM` checkpoints for this fold.

## Training and checkpoint reuse

- `warm-start` and `cold-cm-drug` reuse the three Step 04 Tahoe fold 0 EMA
  checkpoints. `cold-drug` reruns that Tahoe fold and recipe for the three
  seeds after applying the cross-source drug exclusion.
- Within each split and seed, `Tahoe` predicts tenv directly from the selected
  Tahoe checkpoint and `Transfer` initializes from that same checkpoint.
  Rebuild the recipient chemical lookup and dose indices for the tenv input
  and apply the same dose mapping to both arms.
- `Transfer` follows the corresponding `CM` recipe in its tenv stage. Train
  the full model and initialize a fresh optimizer, scheduler and EMA tracker
  from the transferred model parameters. The source and target stages use
  their respective validation data and training recipes.
- `Mix` directly pools the eligible Tahoe and tenv training groups and runs
  one random-init fit. Apply the existing group sampler to the combined pool
  with neither source balancing nor source reweighting. Keep each source's
  group boundaries and matched-control preprocessing. Select the checkpoint
  on the corresponding tenv validation set.

## Evaluation

Within each split regime, the primary contrast pairs `Transfer` and `CM`
by fold and seed. Compute the difference in their equally weighted positive
and negative AP@100 and report the individual paired differences and mean.
Retain all eight Step 04 metrics and both directions of the held-out PR curves.
Seed variation describes repeated fits on a fixed split.
