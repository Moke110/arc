# Separate transcriptome reconstruction from perturbation prediction

**Status:** scope, baseline families and evaluation structure agreed
(2026-11-10). The remaining benchmark settings are under design.

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

The primary task reconstructs each complete measured expression profile
from itself, using the existing cross-attention architecture as the main
model. It includes treated and control profiles, uses the existing
`log1p(CP10K)` expression scale, and uses only genes represented in both
corpora. Any future masked-expression recovery belongs to a separate
experiment.

Use the frozen CM gene order, verifying both sources' native columns by gene
identity before projection. Tenv resolves by HGNC ID. For Tahoe, retain the
exact Ensembl match where it agrees with the identity evidence. If that ID is
absent, use the unique native column linked to the same HGNC ID by the pinned
Step 02 crosswalk. Report any missing, conflicting or ambiguous mappings
before releasing the input. The metadata audit found native columns for all
18,164 CM genes in both sources. FAM153B (HGNC:27323) maps to
`ENSG00000182230` in tenv and `ENSG00000289731` (token `61482`) in Tahoe.
The earlier exact-ID-only projection missed the Tahoe column despite observed
counts. Rebuild the expression release by projecting raw counts before
normalization. All methods then use the measured FAM153B values in both
sources, with observed zero counts retained.

The primary comparison fits each method on the same joint Tahoe and tenv
training membership and reports each source separately. The `mean` reference
is a single vector, with each gene averaged over all joint training profiles.
Source-specific fits are a separate diagnostic.

Freeze disjoint training, validation and test plate memberships independently
within each source. Select epochs and hyperparameters on validation plates
and report the final comparison on test plates. Fit the primary model anew
for this comparison. The historical full-data checkpoint includes the
profiles eligible for the new test split and remains evidence of the earlier
fit.

This ADR records decisions ahead of the directory migration and benchmark.
Historical evidence retains its recorded paths under
[ADR 0030](0030-external-data-provenance-and-step-layout.md).
