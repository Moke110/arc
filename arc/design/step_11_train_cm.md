# Step 11 — Phase 2: Cardiomyocyte Joint Training

## Goal

Fine-tune the full compositional autoencoder (M1–M5) on all 6 cardiomyocyte datasets jointly, with MMD regularization for cross-platform integration and covariate embeddings for context — producing the trained cardiomyocyte drug-perturbation model.

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

- **Gene space:** ~10,000 genes (shared cardiomyocyte gene space, final count after TempO-Seq intersection from step_04). M5 decoder output dimension = final shared gene space size.
- **Covariates for M4:** {cell_type (2), platform (5), study (6), donor (~10+)}.
- **Objective:** Reconstruct cardiomyocyte Δ given control expression + drug/chemical SMILES + dose/concentration + covariates.
- **Loss:** L_recon (Gaussian NLL on Δ) + β_KL · L_KL(M1) + β_MMD · L_MMD(M1, across 5 platforms / 10 pairs).
- **MMD regularization:** computed across platform/study groups within each minibatch. This is the primary cross-platform integration mechanism. With 5 platforms, MMD spans C(5,2)=10 platform pairs.
- **Data:** All 6 cardiomyocyte datasets jointly (train split only — held-out drugs/chemicals excluded). ~6,500+ samples (subject to SMILES filtering of GSE262419 chemicals).
- **Platform-balanced sampling:** Critical — GSE262419's ~4,558 samples are ~70% of total. Without balanced sampling, TempO-Seq dominates every minibatch and MMD is ineffective. Each minibatch should contain approximately equal samples per platform.

### 11.3 Training loop

- Optimizer: Adam or AdamW (specific learning rate, scheduler left open for tuning).
- Batch size: can increase from 32–64 to 64–128 given ~6.5k samples (vs. ~2k previously). Ensure multiple platforms per batch for MMD — platform-balanced sampling is mandatory.
- Epochs: train until convergence or early stopping. With more data, each epoch is slower but convergence may be faster in terms of epochs. Monitor validation loss on a held-out fraction of training samples (not the unseen-drug test set).
- Checkpointing: save model weights periodically to shared checkpoint storage.
- Monitor: training loss, validation loss, MMD loss, per-platform reconstruction loss (5 platforms — detect if TempO-Seq or any platform dominates or underperforms).

### 11.4 Hyperparameter considerations

- **β_MMD (MMD weight):** Critical hyperparameter. Too strong → z_basal loses biological signal. Too weak → platform effects leak. Tune via validation loss. With 5 platforms and 10 MMD pairs, the total MMD loss is a sum over more terms — consider normalizing by the number of pairs.
- **β_KL (KL weight):** Controls VAE regularization strength. Too strong → posterior collapses. Too weak → poor latent structure. Use a small value (e.g., 5e-5 range).
- **M2 freeze vs. fine-tune:** With ~6.5k samples (vs. ~2k previously), fine-tuning M2 is less risky than before. The larger dataset better supports fine-tuning the drug encoder. Reconsider the freeze decision — monitor validation loss with M2 frozen vs. fine-tuned. Default: fine-tune M2 with a small learning rate, but keep the frozen option as a fallback if overfitting is observed.
- **Learning rate scheduling:** Consider warmup + cosine decay or ReduceLROnPlateau.
- **Platform-balanced sampling rate:** The degree of balancing affects training dynamics. Too aggressive (exactly equal per platform) wastes TempO-Seq samples; too lax allows TempO-Seq dominance. Consider a weighted sampler that approximately balances but doesn't discard excessive data.

### 11.5 Training monitoring

Track during training:
- Total loss, reconstruction loss, KL loss, MMD loss (all on training and validation)
- Per-platform reconstruction loss (5 platforms — detect platform imbalance, especially TempO-Seq dominance)
- Per-study reconstruction loss (6 studies — detect study-specific issues)
- Held-out drug/chemical prediction quality (monitor, not optimize — the unseen-drug test set is NOT used for model selection)
- Gradient norms per module (detect vanishing/exploding gradients)
- MMD loss per platform pair (10 pairs — identify which platform pairs are hardest to align)

## Unit tests

1. **Training loop stability:** Verify the training loop runs for at least 50 epochs without NaN/Inf. Verify loss decreases over the first 10 epochs.
2. **MMD active:** Verify MMD loss is computed across all 5 platform groups (10 pairs) in each minibatch. Verify MMD is non-zero when multiple platforms are present.
3. **Held-out drug/chemical exclusion:** Verify no held-out drug/chemical (from step_05 split) appears in the training data. This is a hard constraint.
4. **Checkpoint integrity:** Verify saved checkpoint contains all M1–M5 weights. Verify checkpoint is loadable from a fresh process.
5. **Per-platform loss tracking:** Verify per-platform reconstruction loss is logged for all 5 platforms. Verify no platform has consistently zero loss (would indicate it's not in the training data).
6. **M2 freeze behavior:** If M2 is frozen, verify M2 gradients are zero and M2 weights don't change during training. If fine-tuned, verify M2 gradients are non-zero.
7. **Output dimensionality:** Verify M5 decoder output matches the cardiomyocyte gene space (final shared gene count from step_04, not necessarily exactly 10,000).
8. **Platform-balanced sampling:** Verify that minibatches contain samples from multiple platforms (not dominated by TempO-Seq). Verify no platform exceeds ~50% of any minibatch.

## Integration test

Validate the trained model is evaluation-ready:
- The trained model produces predicted_Δ for any cardiomyocyte drug/chemical (including held-out drugs/chemicals) given control expression + SMILES + dose/concentration + covariates.
- The trained model's predicted_treated = control_expression + predicted_Δ has the correct shape and value range.
- Per-platform reconstruction loss is balanced across 5 platforms (no platform has loss >2x the average — would indicate poor cross-platform integration).
- The model produces non-trivial predictions for held-out drugs/chemicals (predicted_Δ is not all zeros — would indicate the model learned to ignore the drug input).
- Training and validation loss curves show convergence (no severe overfitting — validation loss doesn't diverge from training loss).
- The final checkpoint is saved to shared storage and loadable.

## Acceptance criteria

- [ ] All 6 cardiomyocyte datasets used in training (train split only)
- [ ] M2 initialized from Phase 1 (frozen or fine-tuned, decision documented)
- [ ] MMD regularization active across 5 platforms (10 pairs)
- [ ] Platform-balanced sampling implemented (no TempO-Seq dominance)
- [ ] No held-out drugs/chemicals in training data (hard constraint)
- [ ] Training converges (loss decreases, no NaN/Inf)
- [ ] Per-platform loss balanced across 5 platforms (no platform >2x average)
- [ ] Model produces non-trivial predictions for held-out drugs/chemicals
- [ ] Final checkpoint saved and loadable

## Risks & notes

- GSE262419's ~4,558 samples represent ~70% of total training data. Without platform-balanced sampling, the model will overfit to TempO-Seq characteristics and fail to generalize to RNA-seq platforms. Platform-balanced sampling is the single most important training design decision in this step.
- The larger dataset (~6.5k vs. ~2k) reduces overfitting risk and supports fine-tuning M2. However, the increase is primarily from one dataset (GSE262419) on one platform (TempO-Seq) — chemical diversity increases but platform diversity does not.
- GSE262419's concentration-response design means one chemical contributes multiple samples (one per concentration). This inflates sample count but not chemical diversity. The model sees more Δ profiles per chemical, which helps learn concentration-dependent effects, but the unseen-drug generalization is still limited by the number of unique chemicals (~500+ after SMILES filtering).
- Given Ahlmann-Eltze and PerturBench results, there is a real risk the model does not beat linear baselines. This step must be followed by honest baseline comparison (step_12).
- MMD hyperparameter (β_MMD) is the most critical tuning parameter for cross-platform integration. With 10 platform pairs (vs. 6 previously), the total MMD loss is larger — consider normalizing by the number of pairs to keep β_MMD on the same scale.
- The donor embedding table (M4) may overfit with ~10+ donors and ~6.5k samples. Consider L2 regularization on donor embeddings or reducing donor embedding dimension.
- If the model fails to converge or produces trivial predictions (all Δ ≈ 0), this indicates the compositional structure or MMD regularization is too strong. Reduce β_MMD or β_KL and retry.
- Environmental chemicals from GSE262419 may have very different SMILES structures from clinical drugs. The pretrained M2 encoder (trained on L1000 drug-like compounds) may produce poor embeddings for environmental chemicals. Monitor per-class reconstruction loss for environmental chemicals specifically.
