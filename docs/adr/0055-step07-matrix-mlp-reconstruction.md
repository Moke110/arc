# Compare vector and matrix MLP reconstruction in Step 07

**Status:** accepted (2026-11-10).

Compare the existing vector MLP as `MLPv_512` with a new `MLPm_512`.
Both take and reconstruct one measured profile on the frozen CM gene axis.
`MLPv_512` has 512 latent scalars. `MLPm_512` has 512 tokens and 128
channels per token, for 65,536 latent scalars, matching `AAE_512_SA`.
Keeping the same shape as the cross-attention model lets us compare how
the two models mix gene and latent tokens at the same latent capacity.

Keep the gene embedding and expression readout of `AAE_512_SA`.
Use an MLP shared across the 128 channels to map the gene token axis
to 512 latent tokens. Use a separate channel-shared MLP to map the
latent token axis back to the ordered genes. These learned maps use
the same weights for every profile. The AAE uses attention weights
that depend on its input. Record each model's parameter count alongside
its latent shape when interpreting the comparison.

The revised figure retains only `PCA_512` from the PCA sweep and
`MLPv_512` from the vector MLP sweep. Keep `mean`, the three AAE
variants and `MatrixAE_64x8` alongside the new matrix MLP. Use the
existing `MLP_512` artifact for `MLPv_512`. Preserve the historical
artifact identifiers and seals for replay. PCA continues to use its
512 principal-component scores as a vector.

Use the frozen Step 07 split, validation-only checkpoint selection and
source-specific reconstruction metrics. Commit the revised selection seal
before held-out scoring and independently replay the results.

This revises the active comparison roster in
[ADR 0052](0052-stage2-autoencoder-reconstruction-boundary.md) and extends
the comparison in
[ADR 0054](0054-step07-separate-reconstruction-and-nuisance-tokens.md).
