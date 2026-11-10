# Separate expression reconstruction and nuisance tokens in Step 07

**Status:** accepted (2026-11-10).

Use one gene-to-latent cross-attention layer followed by two latent
self-attention layers. Retain the fixed gene positions, 128 channels per
token and eight attention heads from the published cross-attention baseline.
Compare two variants using the existing gene-query decoder design:

- `AAE_512` uses a 512 × 128 matrix latent representation and decodes all
  512 rows. Compared with the published cross-attention baseline, it adds
  the two encoder self-attention layers.
- `AAE_512+16` uses a 528 × 128 matrix latent representation. The decoder
  reads the first 512 rows. The final 16 are exploratory nuisance tokens.

Both variants use the two pre-normalized Transformer blocks from `MatrixAE`
in the encoder, with residual self-attention and a residual GELU feedforward
network. Each block attends across every latent token, including both token
groups in `AAE_512+16`. The extra queries can learn through their effect on
the reconstruction tokens, which the decoder reads.

Both variants reconstruct the measured expression profile on the frozen CM
axis using full-axis MSE as the sole training objective. Use the existing AAE
fit recipe and ADR 0052's frozen split and source-specific scoring. The first
experiment leaves the specific sources of nuisance variation unspecified.
Reconstruction scores measure the effect of the additional tokens on
reconstruction quality. Their biological interpretation requires separate
evidence.

Label the published cross-attention baseline `AAE_512_SA` in the new
comparison. Its encoder uses one cross-attention layer and zero
self-attention layers. Retain its historical `AAE_512` artifact identifier
for sealed replay. Give the new `AAE_512` variant a distinct artifact
identifier so both can appear in one comparison. Compare the baseline with
the new `AAE_512` to measure the effect of adding the two Transformer blocks.
Compare the two new variants to measure the effect of the extra 16 tokens.

## See

- [ADR 0052](0052-stage2-autoencoder-reconstruction-boundary.md)
- [ADR 0053](0053-step07-compact-matrix-autoencoder.md)
- [Step 07](../../stage2/07_autoencoder/README.md)
