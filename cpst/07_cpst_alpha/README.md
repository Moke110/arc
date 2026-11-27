# Step 07 CarPS-Transcriptome alpha

The fixed alpha is the [Step 06.4 shared signed Δ arm](../06_tsg_shape/06.4_round4/README.md#formal-results) selected in [ADR 0099](../../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md). It uses the ordered tenv Denoise-DEG 5k axis, measured same-plate DMSO and frozen KPGT. Separate 100-dimensional Gaussian encoders model measured control and measured signed Δ. Given a sampled control latent, the native Gaussian prediction head takes a 400-value chemical projection and 32 repeats of the train-fold standardized log-dose. The shared signed decoder maps the measured and predicted Δ latents to the same ordered genes.

The forward path returns the control reconstruction, both control Gaussian parameters, treated prediction `X̂2 = X1 + Δ̂`, both predicted Gaussian parameters and the sampled predicted latent. The objective combines control reconstruction MSE, measured Δ reconstruction MSE, prediction Δ MSE, 1.75 times prediction Δ L1 and 0.1 times the sum of both Gaussian prior KL terms and the prediction-to-measured-Δ Gaussian KL term. Both predicted and measured Δ subtract the same measured control.

## Model API

Add `cpst/07_cpst_alpha/output/api` to the Python import path, then import `AlphaAdapter` from `cpst_alpha`. Construct it with the [pinned TranSiGen source](../../pert_model/TranSiGen/author_log.json) and the training fold's dose transform, then call `build(training_seed, device)`. The adapter's `predict(model, inputs, rows)` accepts the matched control and raw dose supported by the [Step 04 input interface](../04_tsg_axis/README.md), packs the standardized scalar after the chemical vector and returns sampled treated expression in float32. Control and measured responses use the same 5,000 ordered genes. The constructor initializes a model. The public verifier loads the accepted weights and checks their predictions.

The KPGT default uses its native 2,304 inputs and original 400-output projection. To supply another width, pass `chemical_width`. Passing `projection_factory` supplies an experiment-owned callable from `(native_width, 400)` to a CPU `torch.nn.Linear`; the standard variant also uses an affine projection. The adapter constructs the native model first, replaces only that projection under isolated RNG state and moves the new layer to the model device. The variant input reader supplies the appropriate native-width vectors in frozen condition order. Keep this configuration and its outputs with the experiment.

## Checkpoint and replay evidence

The [Step 06.4 comparison](../06_tsg_shape/06.4_round4/README.md#formal-results) reports the fixed alpha as `XΔ-XΔ`. Its [held-out scores](../06_tsg_shape/06.4_round4/figure/metrics.png) and [fold–seed differences](../06_tsg_shape/06.4_round4/figure/metrics-paired.png) use the accepted Step 06.4 fit lineage.

The [release manifest](output/release.json) pins the public API, verifier, runtime dependencies, KPGT defaults and original Step 06.4 revision and source digests. The [parity receipt](log/362-alpha.json) and [resource record](log/362-resources.json) record the capped CUDA replay. The parity receipt binds the registration and checkpoint hashes to exact validation predictions and objective at sampling seed 22617 and exact held-out expression and signed Δ predictions at sampling seed 100042. All nine coordinates cover three cold-drug folds and training seeds 17, 42 and 2026.

From the materialized checkout with the `arc` environment active, use the normal verifier. It reads the sealed Step 06.4 checkpoints and receipts in place and compares the replay with the published parity receipt. Use the original Git history for the pinned Step 06.4 code hashes. Live alternate paths for the same frozen inputs can be passed with `--round-four`, `--training`, `--representations` and `--source`.

```bash
free -g
systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
  env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
  python cpst/07_cpst_alpha/execution/verify_alpha.py
```

Follow the [execution limits](../../docs/agents/execution-environment.md) and record GPU, WSL and Windows memory with the [ten-second coordinator sampler](../06_tsg_shape/06.1_round1/execution/sample_resources.py). See the [pinned old-layout procedure](../06_tsg_shape/06.4_round4/README.md#replay-a-sealed-step-064-fit) for standalone replay through the original Step 06.4 command.
