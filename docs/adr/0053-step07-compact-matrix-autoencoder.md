# Compare a compact matrix latent representation in Step 07

**Status:** accepted (2026-11-10).

Add `MatrixAE_64x8` as a Step 07 transcriptome autoencoder. Its deterministic
bottleneck is 64 latent tokens with eight channels per token. This gives it
the same 512-scalar budget as `PCA_512` and `MLP_512`. Keep the existing
`AAE_512` as the 512-token, 65,536-scalar attention reference.

Follow the scLDM autoencoder's learned pooling and gene-query unpooling.
Learn a gene identity embedding and add a projected expression value at
each measured gene. Use 64 learned queries to pool the full gene axis,
including measured zeros. Apply two latent Transformer blocks before
projecting to eight channels per token. Expand the matrix back to the
internal width of 128, apply two decoder latent Transformer blocks and
reconstruct the full axis using the same gene identity embeddings as
queries. The decoder receives expression information solely through the
64 × 8 bottleneck. Use a linear readout and full-axis MSE on the Step 07
`log1p(CP10K)` scale. Diffusion, KL, stochastic sampling and masked
reconstruction are outside this comparison.

Retain ADR 0052's frozen Tahoe and tenv plate memberships, joint shuffled
training and source-specific metrics. Use the AAE's batch size of 32, AdamW
recipe, 10,000-optimizer-update cap and validation every full epoch and at
the cap. Select the earliest minimum tenv validation MSE and stop after
three consecutive checks without improvement. Bind the selected artifact
to its validation evidence in a committed seal before test scoring. Replay
the held-out scores independently and include the new entrant in the
existing metric figure alongside the previous artifacts.

As a Step 07 exception to [ADR 0030](0030-external-data-provenance-and-step-layout.md),
keep the published benchmark figures at `07_autoencoder/output/*.png` and
the fitted artifacts at `07_autoencoder/output/models/`. Each fitted artifact
directory keeps its tracked manifest and, where applicable, epoch history
together with the ignored local weights or fitted arrays named in the manifest.
This preserves the published figure and replay paths. ADR 0030 continues to
govern the other Stage 2 steps.

## See

- [ADR 0052](0052-stage2-autoencoder-reconstruction-boundary.md)
- [Step 07](../../stage2/07_autoencoder/README.md)
- [Matrix-latent autoencoder research](../research/stage2-07-matrix-latent-autoencoders.md)
- [scLDM paper, Section 3.2 and Appendix H](https://arxiv.org/html/2511.02986v2)
