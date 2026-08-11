# Step 11 — Phase 2: Cardiomyocyte Joint Training

## Goal

Fine-tune the full compositional autoencoder (M1–M5) on all 4 cardiomyocyte datasets jointly, with MMD regularization for cross-platform integration and covariate embeddings for context — producing the trained cardiomyocyte drug-perturbation model.

## Dependencies

- **Step 4 (Harmonize):** Provides the harmonized cardiomyocyte training data (expression, Δ, covariates, SMILES).
- **Step 5 (Drug Split):** Provides the unseen-drug train/test split (only train drugs used in this step).
- **Step 6 (Impl M1+M4):** Provides the cell-state encoder and covariate embedder.
- **Step 7 (Impl M2):** Provides the drug molecular encoder.
- **Step 8 (Impl M3):** Provides the perturbation module.
- **Step 9 (Impl M5):** Provides the Δ-decoder.
- **Step 10 (L1000 Pretrain):** Provides pretrained M2, M3, M5 checkpoints for initialization.

## Inputs

- Harmonized cardiomyocyte training data (train split only): control expression, Δ profiles, drug SMILES, dose, covariate labels (cell_type, platform, study, donor)
- Shared gene space (~10,000 genes)
- Pretrained M2, M3, M5 checkpoints from Phase 1
- Unseen-drug split (from step_05) — held-out drugs are NOT used in training

## Outputs

- Trained model checkpoint (all M1–M5 weights)
- Training logs: loss curves, MMD loss curves, per-platform loss curves
- Held-out drug prediction samples (for monitoring during training, not final evaluation)

## Behavior

### 11.1 Initialization

- **M2 (drug encoder):** Initialize from Phase 1 checkpoint. Freeze or fine-tune with a small learning rate (decision left open — monitor validation loss).
- **M3 (perturbation module):** Initialize from Phase 1 checkpoint. Fine-tune with normal learning rate.
- **M5 (decoder):** Initialize hidden layers from Phase 1 checkpoint. Reinitialize output layer (978 → ~10,000 genes). Fine-tune.
- **M1 (cell-state encoder):** Initialize from scratch (no pretraining — Ahlmann-Eltze: atlas pretraining gives minimal benefit).
- **M4 (covariate embedder):** Initialize from scratch with fresh cardiomyocyte-specific categories (cell_type, platform, study, donor). L1000 cell-line embeddings are discarded.

### 11.2 Training setup

- **Gene space:** ~10,000 genes (shared cardiomyocyte gene space). M5 decoder output dimension = ~10,000.
- **Covariates for M4:** {cell_type (2), platform (≤4), study (4), donor (~10+)}.
- **Objective:** Reconstruct cardiomyocyte Δ given control expression + drug SMILES + dose + covariates.
- **Loss:** L_recon (Gaussian NLL on Δ) + β_KL · L_KL(M1) + β_MMD · L_MMD(M1, across platforms/studies).
- **MMD regularization:** computed across platform/study groups within each minibatch. This is the primary cross-platform integration mechanism.
- **Data:** All 4 cardiomyocyte datasets jointly (train split only — held-out drugs excluded).

### 11.3 Training loop

- Optimizer: Adam or AdamW (specific learning rate, scheduler left open for tuning).
- Batch size: small (e.g., 32–64) given ~2k samples. Ensure multiple platforms/studies per batch for MMD.
- Epochs: train until convergence or early stopping. Monitor validation loss on a held-out fraction of training samples (not the unseen-drug test set).
- Checkpointing: save model weights periodically to shared checkpoint storage.
- Monitor: training loss, validation loss, MMD loss, per-platform reconstruction loss (to detect if one platform dominates).

### 11.4 Hyperparameter considerations

- **β_MMD (MMD weight):** Critical hyperparameter. Too strong → z_basal loses biological signal. Too weak → platform effects leak. Tune via validation loss.
- **β_KL (KL weight):** Controls VAE regularization strength. Too strong → posterior collapses. Too weak → poor latent structure. Use a small value (e.g., 5e-5 range).
- **M2 freeze vs. fine-tune:** With ~2k samples, fine-tuning M2 risks overfitting. Monitor validation loss with M2 frozen vs. fine-tuned. Default: freeze M2, fine-tune M3+M5.
- **Learning rate scheduling:** Consider warmup + cosine decay or ReduceLROnPlateau.

### 11.5 Training monitoring

Track during training:
- Total loss, reconstruction loss, KL loss, MMD loss (all on training and validation)
- Per-platform reconstruction loss (detect platform imbalance)
- Per-study reconstruction loss (detect study-specific issues)
- Held-out drug prediction quality (monitor, not optimize — the unseen-drug test set is NOT used for model selection)
- Gradient norms per module (detect vanishing/exploding gradients)

## Unit tests

1. **Training loop stability:** Verify the training loop runs for at least 50 epochs without NaN/Inf. Verify loss decreases over the first 10 epochs.
2. **MMD active:** Verify MMD loss is computed across platform/study groups in each minibatch. Verify MMD is non-zero when multiple platforms are present.
3. **Held-out drug exclusion:** Verify no held-out drug (from step_05 split) appears in the training data. This is a hard constraint.
4. **Checkpoint integrity:** Verify saved checkpoint contains all M1–M5 weights. Verify checkpoint is loadable from a fresh process.
5. **Per-platform loss tracking:** Verify per-platform reconstruction loss is logged. Verify no platform has consistently zero loss (would indicate it's not in the training data).
6. **M2 freeze behavior:** If M2 is frozen, verify M2 gradients are zero and M2 weights don't change during training. If fine-tuned, verify M2 gradients are non-zero.
7. **Output dimensionality:** Verify M5 decoder output matches the cardiomyocyte gene space (~10,000 genes), not the L1000 978 genes.

## Integration test

Validate the trained model is evaluation-ready:
- The trained model produces predicted_Δ for any cardiomyocyte drug (including held-out drugs) given control expression + SMILES + dose + covariates.
- The trained model's predicted_treated = control_expression + predicted_Δ has the correct shape and value range.
- Per-platform reconstruction loss is balanced (no platform has loss >2x the average — would indicate poor cross-platform integration).
- The model produces non-trivial predictions for held-out drugs (predicted_Δ is not all zeros — would indicate the model learned to ignore the drug input).
- Training and validation loss curves show convergence (no severe overfitting — validation loss doesn't diverge from training loss).
- The final checkpoint is saved to shared storage and loadable.

## Acceptance criteria

- [ ] All 4 cardiomyocyte datasets used in training (train split only)
- [ ] M2 initialized from Phase 1 (frozen or fine-tuned, decision documented)
- [ ] MMD regularization active across platforms/studies
- [ ] No held-out drugs in training data (hard constraint)
- [ ] Training converges (loss decreases, no NaN/Inf)
- [ ] Per-platform loss balanced (no platform >2x average)
- [ ] Model produces non-trivial predictions for held-out drugs
- [ ] Final checkpoint saved and loadable

## Risks & notes

- Small cardiomyocyte sample size (~2k) is the fundamental constraint. The model may underfit or fail to generalize. The L1000-pretrained drug encoder transfers chemical-space knowledge, but the cardiomyocyte-specific adaptation relies on ~2k samples.
- Given Ahlmann-Eltze and PerturBench results, there is a real risk the model does not beat linear baselines. This step must be followed by honest baseline comparison (step_12).
- MMD hyperparameter (β_MMD) is the most critical tuning parameter for cross-platform integration. If MMD is too weak, the model may learn platform-specific shortcuts; if too strong, it may lose biological signal.
- Freezing M2 vs. fine-tuning is a key decision. With ~2k samples, fine-tuning M2 (~0.5M params) risks overfitting. Default to freezing M2 and fine-tuning only M3 + M5 + M1 + M4.
- The donor embedding table (M4) may overfit with ~10+ donors and ~2k samples. Consider L2 regularization on donor embeddings or reducing donor embedding dimension.
- If the model fails to converge or produces trivial predictions (all Δ ≈ 0), this indicates the compositional structure or MMD regularization is too strong. Reduce β_MMD or β_KL and retry.
