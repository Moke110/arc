# Build a whole-transcriptome Latent XPert predictor

**Status:** proposed (2026-11-09). The step path, model separation, latent token count, conditioning inputs, and frozen decoder are confirmed. Training data, supervision, and evaluation are under design review.

`stage2/07_lxpert/` separates transcriptome reconstruction from perturbation prediction. The transcriptome autoencoder learns to reconstruct profiles on the ordered 18,164-gene CM axis. The encoder uses 512 latent queries that cross-attend to expression tokens. The decoder uses gene-position queries that cross-attend to the latent tokens and reconstruct one expression value per gene. Fixed positional embeddings keyed to the frozen axis provide the sole gene identity signal, so axis order forms part of the model contract.

Latent XPert encodes the baseline transcriptome and passes its 512 tokens through an XPert-style treatment branch with block order `CA + SA + SA + CA`. The branch has two condition tokens, a shared Morgan fingerprint projected by a jointly trained linear adapter and a dose embedding. The frozen decoder maps the predicted latent state to post-treatment expression. Subtracting the baseline in expression space yields the full-axis perturbome.

## See

- [ADR 0030](0030-external-data-provenance-and-step-layout.md) for the Stage 2 step layout.
- [ADR 0042](0042-tenv-model-benchmark-redesign.md) for the ordered CM gene axis and the shared Morgan representation.
