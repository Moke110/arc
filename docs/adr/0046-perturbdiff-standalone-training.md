# Tune PerturbDiff training on Tahoe and tenv in one step

**Status: accepted (2026-11-07).** The user confirmed the step boundary,
input axis, split policy, tuning procedure and deliverables.

After the tenv comparison selected PerturbDiff for follow-up, the user
prioritized training configuration and convergence before testing transfer.
`stage2/04_pdiff_hyperparameter/` owns two standalone training studies,
one for the Tahoe pseudobulk corpus and one for tenv. Both train from random
initialization on the same ordered 18,164-gene CM axis and select configurations
and checkpoints using their own validation data.

Tenv reuses the frozen cold-drug splits from 03. Tahoe also holds out
parent-connectivity blocks, keeping each drug's doses and cell lines together.
For each dataset, profile validation at 250 and 500 optimizer steps in short
runs and use the measured cost to choose the formal cadence. Every validation
event is eligible to replace the best checkpoint. Keep the validation history
and express early-stopping patience in optimizer steps. The search covers the
validation interval, minimum training steps, patience window and minimum
improvement. Keep the learning rate, scheduler, batch size and model
architecture fixed. Each outer fold selects its configuration using its own
training and validation data, then freezes it before held-out test evaluation.
Tenv reuses test memberships already scored in 03 for a paired comparison.

## Validation and checkpoint selection

Follow the pinned PerturbDiff author's monitor
`validation_validation_loss_epoch` and minimize it. The author computes each
batch loss as the mean of the native loss multiplied by its weights, then
takes the arithmetic mean of the validation batch losses. This replaces
03's group-count-weighted reduction. Evaluate and retain the EMA weights
associated with the selected loss, following the author's averaging callback.
Use this validation loss to compare configurations and select each run's
best checkpoint. The tuned minimum-improvement threshold governs early
stopping separately from checkpoint ranking. Only the validation loader
participates in selection.

The monitor, loss and weight semantics come from the author's method.
The pilot determines cadence and early-stopping settings. Each run retains
its best checkpoint and compact loss history. The eight 03 test metrics
remain the separate held-out reporting contract.

## Formal run matrix

Use three cold-drug folds and seeds 17, 42 and 2026 for each dataset, for nine
formal fits per dataset. Run each deterministic reference once per fold.
Retain all nine validation-selected tenv checkpoints.

## Motivation

In the published 03 curves, PerturbDiff's nine recorded best validation
checkpoints are at 6,000–13,000 steps, and all nine runs stop at 30,000 steps.
The 03 protocol validates and considers checkpoints every 1,000 optimizer steps.
Its minimum of 20,000 steps and patience of ten validation events put the
earliest possible stop at step 30,000. These observations motivate revisiting
validation cadence, checkpoint selection and stopping together.

The new `04` supersedes ADR 0033's numbering gap. It consumes the existing
Tahoe and tenv corpora owned by `02_dataset_corpus/`. The retired Tahoe data
step remains folded into the corpus layer.

## Deliverables

- `output/learning_curve.png` has a Tahoe panel and a tenv panel. Training is
  blue, validation is green, and red marks the validation-selected checkpoint.
- `output/metrics.png` has one bar-chart panel per held-out test metric from
  03. Tahoe PerturbDiff is blue, tenv PerturbDiff is green, and the
  `no_change` and `delta_mean` references are gray. Score each reference on
  the same held-out groups as PerturbDiff for its dataset.
- `output/` retains the validation-selected checkpoint of every formal tenv
  run for downstream use.
