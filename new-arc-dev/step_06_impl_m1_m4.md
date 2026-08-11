# Step 6 — Model Impl: M1 (Cell-State Encoder) + M4 (Covariate Embedder)

## Goal

Implement the cell-state encoder (M1) that maps control expression to a platform-invariant basal-state latent, and the covariate embedder (M4) that encodes technical/biological context into a covariate vector — the two modules that together handle cross-platform heterogeneity integration (G1).

## Dependencies

- **Step 0 (Review):** Provides the evidence basis for MMD regularization (instead of adversarial disentanglement) and explicit covariate embedding.

## Inputs

- `control_expression`: G-dimensional vector of log1p-normalized gene expression from a control (DMSO/vehicle) sample
- Covariate indices: {cell_type, platform, study, donor} as integer category indices

## Outputs

- `z_basal`: basal-state latent vector (L-dimensional), platform-invariant
- `z_cov`: covariate vector (C-dimensional), encoding technical/biological context
- Loss contributions: KL divergence (M1 VAE) and MMD regularization (M1 cross-platform invariance)

## Behavior

### M1 — Cell-state encoder

**Functional objective:** Encode the unperturbed cardiomyocyte transcriptome into a compact basal-state latent that captures the cell's baseline biology, stripped of batch/platform identity.

**Input → Output:** control_expression (G-dim) → z_basal (L-dim)

**Key requirements:**
- Must produce a stochastic latent (VAE-style: output both mean and variance, sample via reparameterization).
- Must be platform-invariant: the same control cardiomyocyte should produce the same basal state regardless of measurement platform. This is enforced by MMD regularization across platform/study groups within each minibatch — NOT by an adversarial classifier (per PerturBench evidence that CPA-noAdv > CPA).
- MMD is computed across all pairs of platforms/studies in the minibatch using a mixture of RBF kernels.
- No atlas/foundation-model pretraining for M1 (Ahlmann-Eltze: atlas pretraining gives minimal benefit). Trained from scratch in Phase 2.

**Loss contributions:**
- KL divergence: D_KL(q(z_basal|x) ‖ N(0, I)), with a small weight (β_KL).
- MMD regularization: MMD²(z_basal^platform_A, z_basal^platform_B) across platform/study pairs, with weight β_MMD.

**Contribution to goals:**
- G1 (heterogeneity): MMD regularization is the primary cross-platform integration mechanism. Forces z_basal to discard technical signal and retain biological signal.
- G2 (unseen-drug): A clean, disentangled basal state means the perturbation module only needs to learn drug-specific effects, not re-encode cell identity.

### M4 — Covariate embedder

**Functional objective:** Encode the technical and biological context of each sample into a covariate vector that captures residual platform/batch/donor effects not removed by M1's MMD regularization, so the decoder can account for them without contaminating the basal state.

**Input → Output:** {cell_type, platform, study, donor} category indices → z_cov (C-dim)

**Key requirements:**
- Separate learnable embedding lookup for each covariate (cell_type, platform, study, donor).
- Concatenate or combine all four embeddings into a single covariate vector.
- No standalone loss — trained end-to-end through the decoder's reconstruction loss.
- The MMD loss on M1 creates soft pressure for z_cov to absorb technical variation (since z_basal is forced to be platform-invariant, platform-specific variance must go somewhere — it goes to z_cov).

**Contribution to goals:**
- G1 (heterogeneity): M4 is the explicit context-embedding module that Wei 2025 identified as the key gap. By giving each platform/study/donor its own embedding, the model accounts for systematic differences without forcing the basal state or perturbation vector to encode them.
- G2 (unseen-drug): Minimal direct contribution, but by cleanly separating context from drug effect, it prevents conflating platform artifacts with drug responses.

**Pretraining note:** During Phase 1 (L1000), M4 is trained with L1000 cell-line covariates. These embeddings are discarded in Phase 2 — only the M4 architecture carries over. Phase 2 initializes fresh embeddings for cardiomyocyte-specific categories.

## Unit tests

1. **M1 forward shape:** Given a control_expression input of shape (batch, G), verify z_basal output has shape (batch, L). Verify both mean and variance are output.
2. **M1 reparameterization:** Verify that sampling is differentiable (gradient flows through the reparameterization trick). Verify that with variance → 0, the sampled z_basal → mean.
3. **M1 MMD computation:** Given two groups of z_basal vectors (simulating two platforms), verify MMD is non-negative, zero when the groups are identical, and positive when they differ. Verify MMD uses a mixture of RBF kernels.
4. **M4 forward shape:** Given covariate indices for a batch, verify z_cov output has shape (batch, C).
5. **M4 embedding independence:** Verify that changing one covariate (e.g., platform) while holding others constant changes z_cov. Verify that each covariate has its own embedding table.
6. **Composition:** Verify that z = z_basal + z_pert + z_cov is computed element-wise (additive) and has the correct dimensionality for the decoder input.

## Integration test

Validate M1 + M4 work together for cross-platform integration:
- Given control expression from two simulated "platforms" (same biological state, different technical noise), verify that after MMD regularization training, z_basal distributions across platforms are more similar than before (MMD decreases).
- Given covariate labels for different platforms, verify z_cov captures platform-specific variation (z_cov differs across platforms for the same biological state).
- Verify the combined z = z_basal + z_cov has the correct shape for decoder input.
- Verify that M1 + M4 can process a minibatch with samples from multiple platforms/studies simultaneously.

## Acceptance criteria

- [ ] M1 produces stochastic z_basal with reparameterization
- [ ] M1 MMD regularization computed across platform/study pairs
- [ ] M4 produces z_cov from 4 covariate inputs
- [ ] M4 has separate embedding tables per covariate
- [ ] Additive composition z = z_basal + z_pert + z_cov verified
- [ ] No adversarial classifier in M1 (per PerturBench evidence)

## Risks & notes

- MMD regularization strength (β_MMD) is a key hyperparameter — too strong and z_basal loses biological signal; too weak and platform effects leak through. This is tuned during training (step_11).
- The covariate embedder's embedding dimensions are a design decision left open — they should be large enough to capture platform differences but small enough to avoid overfitting.
- With ~2k cardiomyocyte samples and ~10+ donors, the donor embedding table may overfit. Consider regularization on donor embeddings.
- The decision to use MMD instead of adversarial disentanglement is evidence-based (PerturBench: CPA-noAdv > CPA) but MMD has its own hyperparameters (kernel bandwidths) that need tuning.
