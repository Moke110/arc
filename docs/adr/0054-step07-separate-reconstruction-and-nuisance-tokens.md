# Separate expression reconstruction and nuisance tokens in Step 07

**Status:** proposed (2026-11-10).

Use one gene-to-latent cross-attention layer followed by two latent
self-attention layers. Retain the fixed gene positions, 128 channels per
token and eight attention heads from `AAE_512`. The encoder returns a
528 × 128 matrix latent representation. The existing gene-query decoder
design reads only the first 512 rows. The last 16 rows are nuisance tokens,
designated to capture variation outside the biological signal of interest.

The decoder slice fixes which rows the decoder reads. If the encoder allows
the first 512 rows to attend to the last 16 rows, nuisance token states can
still affect expression reconstruction through encoder self-attention. The
definition of nuisance variation, the reconstruction target, attention
connectivity and training signals are open design decisions. The current
measured-profile reconstruction benchmark is defined by ADR 0052.

## See

- [ADR 0052](0052-stage2-autoencoder-reconstruction-boundary.md)
- [ADR 0053](0053-step07-compact-matrix-autoencoder.md)
- [Step 07](../../stage2/07_autoencoder/README.md)
