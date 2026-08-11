# Step 9 — Model Impl: M5 (Δ-Decoder)

## Goal

Implement the Δ-decoder (M5) that maps the compositional latent vector (basal state + perturbation + covariate) back to a G-dimensional predicted transcriptomic change (Δ), which is added to control expression to produce the predicted post-treatment profile.

## Dependencies

- **Step 0 (Review):** Provides the evidence basis for predicting Δ (not absolute expression) as the platform-invariance mechanism (AetherCell, PRnet).

## Inputs

- `z = z_basal + z_pert + z_cov`: compositional latent vector (L+C-dimensional, additive composition)

## Outputs

- `predicted_Δ`: G-dimensional vector of predicted log-expression changes
- `predicted_treated = control_expression + predicted_Δ`: G-dimensional predicted post-treatment expression
- Loss contribution: Gaussian negative log-likelihood (reconstruction loss on Δ)

## Behavior

### M5 — Δ-decoder

**Functional objective:** Map the compositional latent vector back to a G-dimensional predicted transcriptomic change (Δ), which is added to the control expression to produce the predicted post-treatment profile.

**Input → Output:** z (compositional latent) → predicted_Δ (G-dim)

**Key requirements:**
- Must predict **Δ** (change-from-control), not absolute expression. Final output: predicted_treated = control_expression + predicted_Δ. This is the platform-invariance mechanism — Δ cancels platform-specific systematic effects that dominate absolute expression.
- Must output both a mean (μ_Δ) and a variance (log σ²_Δ) for a Gaussian likelihood, allowing the model to express uncertainty — important for unseen drugs where the model should be less confident.
- The Gaussian NLL reconstruction loss is the primary training objective for the entire model. It flows back through z into M1 (via z_basal), M3 (via z_pert), M4 (via z_cov), and M2 (via z_pert → drug_emb).
- The specific decoder architecture (MLP, residual network, etc.) is left open for further decision. The key constraint is that it must map the compositional latent to a G-dimensional Δ vector with Gaussian likelihood.

**Loss function (primary training objective for the whole model):**
- Gaussian NLL: L_recon = Σ_g [ (μ_Δ,g − Δ_true,g)² / (2σ²_g) + ½ log σ²_g ]
- where Δ_true = log1p(treated) − log1p(control)
- Total loss = L_recon + β_KL · L_KL(M1) + β_MMD · L_MMD(M1)

**Pretraining:** M5 is pretrained in Phase 1 on L1000 (predicting L1000 Δ-profiles on 978 landmark genes). In Phase 2, the decoder output dimension is expanded from 978 to G (~10,000) and fine-tuned on cardiomyocyte data.

**Contribution to goals:**
- G1 (heterogeneity): Predicting in Δ-space (not absolute expression) is the platform-invariance mechanism. Absolute expression is dominated by platform-specific gene-length bias, normalization artifacts, and detection limits. Δ (treated − control within the same experiment) cancels most platform-specific systematic effects, making the prediction target more comparable across platforms.
- G2 (unseen-drug): The decoder must produce a biologically correct Δ for an unseen drug given its perturbation vector. The Gaussian variance head allows the model to express higher uncertainty for unseen drugs.

## Unit tests

1. **Forward shape:** Given z of shape (batch, L+C), verify predicted_Δ has shape (batch, G). Verify both μ_Δ and log σ²_Δ are output.
2. **Gaussian NLL computation:** Given μ_Δ, log σ²_Δ, and Δ_true, verify the NLL is non-negative. Verify that NLL is minimized when μ_Δ = Δ_true and σ²_Δ is small. Verify that NLL increases when μ_Δ deviates from Δ_true.
3. **Δ prediction (not absolute):** Verify that predicted_treated = control_expression + predicted_Δ (element-wise addition). Verify the decoder outputs Δ, not absolute expression.
4. **Variance head:** Verify that log σ²_Δ is a valid log-variance (can be any real number; σ² = exp(log σ²) > 0). Verify that the variance head allows different uncertainty per gene.
5. **Differentiability:** Verify gradients flow from L_recon back through the decoder to z (and thus to M1, M3, M4, M2).
6. **Output dimensionality expansion:** Verify the decoder can be configured for 978 genes (L1000 pretrain) and reconfigured for ~10,000 genes (cardiomyocyte finetune) without architectural changes to the input side.

## Integration test

Validate M5 works in the full model composition:
- Given control_expression → M1 → z_basal, drug_SMILES → M2 → M3 → z_pert, covariates → M4 → z_cov, verify z = z_basal + z_pert + z_cov → M5 → predicted_Δ has shape (batch, G).
- Verify predicted_treated = control_expression + predicted_Δ.
- Verify the total loss (L_recon + β_KL · L_KL + β_MMD · L_MMD) is a scalar that can be backpropagated.
- Verify that for a no-perturbation input (dose=0, z_pert ≈ 0), predicted_Δ is close to zero (the model predicts no change when there is no drug).
- Verify that the decoder can switch between 978-gene output (L1000) and ~10,000-gene output (cardiomyocyte) by changing only the output layer configuration.

## Acceptance criteria

- [ ] M5 maps compositional latent to predicted_Δ (G-dim) with Gaussian likelihood
- [ ] predicted_treated = control_expression + predicted_Δ verified
- [ ] Gaussian NLL loss computed correctly
- [ ] Variance head outputs per-gene uncertainty
- [ ] Gradients flow from L_recon to all upstream modules
- [ ] Output dimensionality configurable (978 for L1000, ~10,000 for cardiomyocyte)

## Risks & notes

- The Δ-space prediction assumes that platform-specific effects cancel in the treated − control difference. This is approximately true for within-experiment comparisons but may not fully hold across platforms with very different gene coverage or dynamic ranges.
- The Gaussian likelihood assumes per-gene independent noise, which is a simplification — gene expression changes are correlated. This is a standard assumption in CPA/chemCPA and is acceptable for a first model.
- The variance head is important for honest uncertainty quantification on unseen drugs, but it can also be exploited by the model to reduce the loss by inflating variance (a form of underfitting). Monitor the average predicted variance during training.
- The output dimensionality expansion (978 → ~10,000) between Phase 1 and Phase 2 means the decoder's output layer is reinitialized. The hidden layers can carry over, but the output projection must be retrained.
