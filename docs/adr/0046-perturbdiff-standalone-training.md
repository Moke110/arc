# Tune PerturbDiff training in separate Tahoe and tenv studies

**Status: proposed (2026-11-07).** The scope and step names are agreed. The
input, split, training and evaluation protocols remain under design.

After the tenv comparison selected PerturbDiff for follow-up, the user
prioritized training configuration and convergence before testing transfer.
`stage2/04_pdiff_tahoe/` and `stage2/05_pdiff_tenv/` will each fit PerturbDiff
from random initialization and use their own source corpus for training,
validation and held-out testing. Each study selects its own checkpoint from
validation. Both studies will validate more frequently than the tenv benchmark,
so early improvements and subsequent plateaus can be measured more closely.

In the published 03 curves, PerturbDiff's nine recorded best validation
checkpoints are at 6,000–13,000 steps, and all nine runs stop at 30,000 steps.
The 03 protocol validates and considers checkpoints every 1,000 optimizer steps.
Its minimum of 20,000 steps and patience of ten validation events put the
earliest possible stop at step 30,000. These observations motivate revisiting
validation cadence, checkpoint selection and stopping together.

The new `04` supersedes ADR 0033's numbering gap. It consumes the existing
Tahoe corpus from `02_dataset_corpus/`. The retired Tahoe data step remains
folded into the corpus layer.
