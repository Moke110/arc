# Separate expression reconstruction and nuisance tokens in Step 07

**Status:** proposed (2026-11-10).

Use one gene-to-latent cross-attention layer followed by two latent
self-attention layers. Retain the fixed gene positions, 128 channels per
token and eight attention heads from the published cross-attention baseline.
Compare two variants using the existing gene-query decoder design:

- `AAE_512` uses a 512 × 128 matrix latent representation and decodes all
  512 rows. Compared with the published cross-attention baseline, it adds
  the two encoder self-attention layers.
- `AAE_512+16` uses a 528 × 128 matrix latent representation. The decoder
  reads the first 512 rows. The final 16 are exploratory nuisance tokens.

Both variants reconstruct the measured expression profile on the frozen CM
axis using full-axis MSE as in ADR 0052. The first experiment leaves the
specific sources of nuisance variation unspecified. Reconstruction scores
measure the effect of the additional tokens on reconstruction quality.

If the encoder allows the first 512 rows to attend to the last 16 rows,
the nuisance token states can affect reconstruction through encoder
self-attention. Attention connectivity, the self-attention block design,
and how to distinguish the new `AAE_512` from the published artifact of the
same name remain open.

## See

- [ADR 0052](0052-stage2-autoencoder-reconstruction-boundary.md)
- [ADR 0053](0053-step07-compact-matrix-autoencoder.md)
- [Step 07](../../stage2/07_autoencoder/README.md)
