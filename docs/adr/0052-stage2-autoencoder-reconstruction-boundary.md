# Separate transcriptome reconstruction from perturbation prediction

**Status:** scope and baseline families agreed (2026-11-10). The benchmark
protocol is under design.

Step 07 becomes `stage2/07_autoencoder/` and owns transcriptome reconstruction
on the frozen CM gene axis. Step 08 owns perturbation prediction and its
integration checks. The Step 08 directory name is deferred. This revises the
step ownership in [ADR 0051](0051-latent-xpert.md), which records the original
implementation and evidence.

Step 07 compares reconstruction with the training-gene mean (`mean`), PCA
with inverse transformation, and a fully connected MLP autoencoder. The
training-gene mean predicts an expression vector estimated per gene from
training profiles alone. The previous figure's matched DMSO and matched DMSO
plus mean drug effect use a different input and belong to the perturbation
prediction evaluation. The rationale and literature are in the
[reconstruction research note](../research/stage2-07-transcriptome-reconstruction-baselines.md).

The benchmark design still needs decisions on plate membership, training by
source, and the role of the existing cross-attention autoencoder. This ADR
records the agreed scope ahead of the directory migration and benchmark.
Historical evidence retains its recorded paths under
[ADR 0030](0030-external-data-provenance-and-step-layout.md).
