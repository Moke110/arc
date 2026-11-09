# Build a whole-transcriptome Latent XPert predictor

**Status:** proposed (2026-11-09). The gene representation, latent token count, and chemical representation are confirmed. Step numbering, the prediction and training contract, and evaluation scope are under design review.

Latent XPert uses the ordered 18,164-gene CM axis. Each gene's expression token receives a fixed positional embedding keyed to its index in that axis. The axis order therefore forms part of the model's input and output contract. The first version uses position as the sole gene identity signal and takes its chemical input from the shared Morgan fingerprint with a jointly trained linear adapter.

The proposed architecture uses 256 latent cell tokens. Learned queries read the full gene axis through cross-attention and provide a common input to an untreated branch and a drug-conditioned branch. Gene-position queries read the latent representations to produce a full-axis perturbome. This confines dense cell self-attention to the latent sequence while preserving a prediction for each gene in the ordered axis.

## See

- [ADR 0030](0030-external-data-provenance-and-step-layout.md) for the Stage 2 step layout.
- [ADR 0042](0042-tenv-model-benchmark-redesign.md) for the ordered CM gene axis and the shared Morgan representation.
- [ADR 0050](0050-tahoe-transfer-across-stage2-models.md) for the separate proposed multi-model transfer study at `stage2/07_models_transfer/`.
