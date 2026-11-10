# Separate transcriptome reconstruction from perturbation prediction

**Status:** accepted (2026-11-10).

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
normalization. For each profile, divide each projected count by the sum
over the 18,164 shared genes, multiply by 10,000, then apply `log1p`.
This gives both sources the same denominator gene set and removes a
per-profile multiplicative count scale. All methods use the measured
FAM153B values in both sources and retain observed zero counts.

The primary comparison fits each method on the same joint Tahoe and tenv
training membership and reports each source separately. The `mean` reference
is a single vector, with each gene averaged over all joint training profiles.
Source-specific fits are a separate diagnostic. Compare PCA and MLP AE at
32, 128 and 512 latent scalars. Report the main model's 512 latent tokens
with 128 channels each separately so token count is distinguishable from
latent scalar count. [ADR 0055](0055-step07-matrix-mlp-reconstruction.md)
revises the active comparison roster while preserving the sealed sweep.

Freeze disjoint training, validation and test plate memberships independently
within each source. Keep the existing validation plates, then use the next
equally sized block of plates from the existing source-specific hash ordering
with seed 20261109 as test. The rest form the training split.

| Source | Train plates / profiles | Validation plates / profiles | Test plates / profiles |
| --- | ---: | ---: | ---: |
| Tahoe | 10 / 47,822 | 2 / 9,600 | 2 / 9,596 |
| tenv | 10 / 2,810 | 3 / 843 | 3 / 825 |

Fit all methods anew on the training split. For epoch-based models, visit
each training profile once per epoch at the natural source ratio. Select
epochs and configurations by tenv validation MSE, and report Tahoe validation
scores separately. Report test PCC and RMSE per profile and average them
equally within each source. Report squared reconstruction error relative to
the training-gene mean on the same test profiles to quantify how much
profile-specific variation is preserved.

Retain the exact fitted parameters used for test scoring as local artifacts.
Publish their manifests and digests with the split, gene axis, and scores in
Git. Sealed replay requires access to the exact local payloads identified by
those digests. The historical full-data checkpoint remains local evidence of
the earlier fit. Any future Step 08 refit must follow the test exclusions of its
perturbation evaluation. Historical evidence retains its recorded paths under
[ADR 0030](0030-external-data-provenance-and-step-layout.md).
