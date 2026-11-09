# Build a whole-transcriptome Latent XPert predictor

**Status:** accepted (2026-11-09).

The [original implementation and evidence](../../stage2/_archive/07_lxpert/README.md)
are described below.
[ADR 0052](0052-stage2-autoencoder-reconstruction-boundary.md) records the revised
step boundary and reconstruction references.

`stage2/07_lxpert/` separates transcriptome reconstruction from perturbation prediction. The transcriptome autoencoder learns to reconstruct profiles on the ordered 18,164-gene CM axis. The encoder uses 512 latent queries that cross-attend to expression tokens. The decoder uses gene-position queries that cross-attend to the latent tokens and reconstruct one expression value per gene. Fixed positional embeddings keyed to the frozen axis provide the sole gene identity signal, so axis order forms part of the model contract.

Use 128 channels per latent token and eight attention heads. The bottleneck limits the attention sequence to 512 tokens. It contains 65,536 scalar values per profile.

The autoencoder training population combines the Step 02 Tahoe pseudobulk corpus and all tenv treated and control transcriptomes. Technical blanks do not contain a transcriptome and are excluded. Project raw counts onto the frozen CM axis before calculating CP10K and `log1p`. The reconstruction task uses expression profiles alone and does not apply the perturbation predictor's chemical eligibility or matched-control requirements. This population includes the treated profiles reserved for the earlier tenv cold-drug benchmark. The resulting autoencoder weights are for reconstruction and integration testing. A future comparable cold-drug predictor will require autoencoder training that excludes its held-out treated profiles.

During development, reserve validation plates independently within Tahoe and tenv. Each training epoch visits every eligible training profile once at the natural source ratio. Select the epoch count by tenv validation reconstruction MSE and report Tahoe validation MSE separately. Refit from random initialization at the selected epoch count using every eligible profile and publish that checkpoint. The plate-held-out scores describe the development fit.

Latent XPert encodes the baseline transcriptome and passes its 512 tokens through an XPert-style treatment branch with block order `CA + SA + SA + CA`. The branch has two condition tokens, a shared Morgan fingerprint projected by a jointly trained linear adapter and an affine projection of `log10(d / 0.05 µM)` for measured positive dose. The frozen decoder maps the predicted latent state to post-treatment expression. Subtracting the measured same-plate DMSO baseline yields the full-axis perturbome. The downstream training contract freezes both encoder and decoder and lets the full-axis Δ MSE and Pearson loss backpropagate through the decoder into the predictor and its condition adapters.

Step 07 formally trains and releases the transcriptome autoencoder. The perturbation predictor receives only an integration smoke test in this step. Formal perturbation training, held-out perturbation scoring, and comparison against PerturbDiff are outside this step's scope.

Validate the data release, plate split, bounded fit and checkpoint reload, and full-axis predictor integration through the Step 07 command-line interface.

## See

- [ADR 0030](0030-external-data-provenance-and-step-layout.md) for the Stage 2 step layout.
- [ADR 0033](0033-stage2-dataset-corpus.md) for the Tahoe and tenv source corpora.
- [ADR 0042](0042-tenv-model-benchmark-redesign.md) for the ordered CM gene axis and the shared Morgan representation.
- [ADR 0049](0049-perturbdiff-transfer-optimization-comparisons.md) for continuous log-dose conditioning.
