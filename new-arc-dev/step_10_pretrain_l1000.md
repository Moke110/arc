# Step 10 — Phase 1: L1000 Drug-Encoder Pretrain

## Goal

Pretrain the drug molecular encoder (M2), perturbation module (M3), and Δ-decoder (M5) on the LINCS L1000 corpus (~1.3M drug-perturbation profiles, ~20k drugs, ~70 cell lines) to give the drug encoder broad chemical-space coverage so unseen drugs can be encoded meaningfully.

## Dependencies

- **Step 3 (Retrieve L1000):** Provides the L1000 Δ-profile matrix, metadata, and SMILES mappings.
- **Step 4 (Harmonize):** Provides the harmonized L1000 data on 978 landmark genes with covariate labels.
- **Step 7 (Impl M2):** Provides the drug molecular encoder.
- **Step 8 (Impl M3):** Provides the perturbation module.
- **Step 9 (Impl M5):** Provides the Δ-decoder (configured for 978-gene output).

## Inputs

- L1000 Δ-profile matrix (~1.3M samples × 978 genes)
- L1000 metadata (drug SMILES, dose, timepoint, cell line, batch)
- L1000 cell-line covariate labels (for M4, trained but discarded in Phase 2)
- Initialized model modules: M2, M3, M5 (and M1, M4 trained with L1000 covariates but weights discarded)

## Outputs

- Pretrained M2 checkpoint (drug molecular encoder weights)
- Pretrained M3 checkpoint (perturbation module weights)
- Pretrained M5 checkpoint (Δ-decoder weights, 978-gene output layer)
- Training logs: loss curves, embedding quality metrics
- Drug embedding quality report: structure-similarity preservation on held-out L1000 drugs

## Behavior

### 10.1 Training setup

- **Model:** Full compositional autoencoder (M1–M5), but only M2, M3, M5 weights are saved for Phase 2. M1 and M4 are trained with L1000 covariates but their weights are discarded.
- **Gene space:** 978 L1000 landmark genes. M5 decoder output dimension = 978.
- **Covariates for M4:** cell_line (~70 categories), batch. These embeddings are discarded in Phase 2.
- **Objective:** Reconstruct L1000 Δ given drug SMILES + dose + cell-line covariate.
- **Loss:** L_recon (Gaussian NLL on Δ) + β_KL · L_KL(M1) + β_MMD · L_MMD(M1, across cell lines).

### 10.2 Training loop

- Optimizer: Adam or AdamW (specific learning rate, scheduler left open for tuning).
- Batch size: large (e.g., 256–512) to ensure multiple cell lines per batch for MMD computation.
- Epochs: train until loss converges or early stopping triggers (specific epoch count left open).
- MMD regularization: computed across cell-line groups within each minibatch.
- Checkpointing: save model weights periodically to shared checkpoint storage.
- Validation: hold out a fraction of L1000 drugs (unseen-drug split within L1000) to monitor embedding quality during pretraining.

### 10.3 Embedding quality assessment

After pretraining, evaluate the drug encoder's embedding quality:
- **Structure-similarity preservation:** For held-out L1000 drugs, verify that structurally similar drugs produce similar embeddings (cosine similarity).
- **Drug recovery:** For a held-out drug, verify that its nearest neighbors in embedding space are structurally similar drugs (Tanimoto similarity of Morgan fingerprints).
- **Δ reconstruction:** Verify the pretrained model can reconstruct L1000 Δ-profiles for held-out drugs with non-trivial accuracy (better than Δ=0 baseline).

### 10.4 Checkpoint saving

Save the following to shared checkpoint storage for Phase 2 initialization:
- M2 weights (drug molecular encoder) — the most important transfer
- M3 weights (perturbation module)
- M5 weights (decoder hidden layers; output layer will be reinitialized for cardiomyocyte gene space)
- Training logs and embedding quality report

## Unit tests

1. **Training loop stability:** Verify the training loop runs for at least 100 batches without NaN/Inf in loss or gradients. Verify loss decreases over the first 100 batches.
2. **MMD active during pretraining:** Verify MMD loss is computed across cell-line groups in each minibatch. Verify MMD is non-zero when multiple cell lines are present in the batch.
3. **Checkpoint integrity:** Verify saved checkpoint files are loadable and contain the expected parameter keys (M2, M3, M5 weights). Verify checkpoint file sizes are non-zero.
4. **Loss convergence:** Verify the training loss decreases over epochs (not flat or increasing). Verify validation loss is tracked separately.
5. **Batch composition:** Verify each minibatch contains samples from multiple cell lines (needed for MMD). Flag if any batch has only one cell line.
6. **Gradient flow:** Verify gradients are computed for M2, M3, M5 parameters after each batch. Verify no parameter has consistently zero gradients (would indicate a broken computational graph).

## Integration test

Validate the pretrained encoder is Phase-2-ready:
- M2 checkpoint loads correctly and produces embeddings for cardiomyocyte drugs (drugs not in L1000 should still produce valid embeddings if their SMILES are valid).
- M3 checkpoint loads correctly and accepts M2 embeddings + dose.
- M5 checkpoint loads correctly and produces 978-dimensional Δ predictions.
- The embedding quality report shows structure-similarity preservation (held-out L1000 drugs cluster by structural similarity in embedding space).
- The pretrained model reconstructs L1000 Δ-profiles for held-out drugs better than the Δ=0 baseline.
- All checkpoints are saved to shared storage and are loadable from a fresh process.

## Acceptance criteria

- [ ] Training loop runs to completion (or early stopping) without NaN/Inf
- [ ] Loss converges (training and validation loss decrease)
- [ ] MMD regularization active across cell-line groups
- [ ] M2, M3, M5 checkpoints saved to shared storage
- [ ] Embedding quality report shows structure-similarity preservation
- [ ] Pretrained model beats Δ=0 baseline on held-out L1000 drugs
- [ ] Checkpoints loadable from a fresh process

## Risks & notes

- L1000 pretraining is the longest single compute step (~4–8h GPU). Budget accordingly and checkpoint frequently.
- L1000 cell lines are all cancer cell lines — drug mechanisms may differ from cardiomyocytes. The pretrain gives chemical-structure-aware embeddings, not cardiomyocyte-specific drug responses. Phase 2 fine-tuning (step_11) adapts the mapping.
- If GPU memory is constrained, reduce batch size or use gradient accumulation. MMD requires multiple cell lines per batch, so batch size cannot be too small.
- The decision to discard M1 and M4 weights from Phase 1 is deliberate — L1000 cell-line covariates are irrelevant to cardiomyocyte training. Only the drug-encoder pathway (M2→M3) and decoder hidden layers (M5) carry useful transfer.
- If L1000 download is infeasible (step_03 risk), this step can be run with the reduced 875-drug signature set, with the tradeoff of less chemical-space coverage and weaker unseen-drug generalization.
