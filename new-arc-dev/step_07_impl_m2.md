# Step 7 — Model Impl: M2 (Drug Molecular Encoder)

## Goal

Implement the drug molecular encoder (M2) that maps a drug's molecular structure (SMILES) into a fixed-dimensional embedding capturing pharmacophore and mechanism-of-action features — the critical module for unseen-drug generalization (G2).

## Dependencies

- **Step 0 (Review):** Provides the evidence basis for molecular-structure-based drug encoding and L1000 pretraining (Ahlmann-Eltze: perturbation-data pretraining is the only transfer that helps).

## Inputs

- `drug_SMILES`: canonical SMILES string for a drug compound

## Outputs

- `drug_emb`: D-dimensional drug molecular embedding

## Behavior

### M2 — Drug molecular encoder

**Functional objective:** Map a drug's molecular structure into a fixed-dimensional embedding that captures its pharmacophore and mechanism-of-action features, enabling the model to encode drugs it has never seen during cardiomyocyte training.

**Input → Output:** drug_SMILES (string) → drug_emb (D-dim)

**Key requirements:**
- Must accept a SMILES string and convert it to a molecular representation suitable for neural network processing.
- Must be structure-aware: drugs with similar molecular structures should produce similar embeddings. This is the mechanism that enables unseen-drug generalization — a held-out drug structurally similar to a training drug will produce a similar embedding and thus a similar predicted perturbation.
- Must be differentiable end-to-end (gradients flow from the reconstruction loss through M3 into M2).
- No standalone loss — trained end-to-end through M3's reconstruction loss.
- The specific graph/network architecture (e.g., message-passing, attention, fingerprint-based) is left open for further decision. The key constraint is that it must preserve structural similarity in the embedding space.

**Pretraining (Phase 1):** M2 is pretrained on LINCS L1000 (~20,000 drugs, ~1.3M profiles). The L1000 pretraining teaches the encoder to map molecular structure → drug-response-relevant features across a vast chemical space. After Phase 1, M2's weights are frozen or fine-tuned with a small learning rate in Phase 2.

**Contribution to goals:**
- G2 (unseen-drug): This is the critical module for unseen-drug generalization. By encoding drugs from molecular structure (not learned drug-ID embeddings), the model can encode any drug with a known SMILES — including held-out drugs never tested on cardiomyocytes. Structure-similarity preservation is the core requirement.
- G1 (heterogeneity): Indirect — the drug embedding is platform-agnostic (molecular structure doesn't depend on assay), providing a stable anchor across platforms.

## Unit tests

1. **SMILES → embedding shape:** Given a valid SMILES string, verify drug_emb output has shape (D,). Given a batch of B SMILES, verify output has shape (B, D).
2. **Structure-similarity preservation:** Given two structurally similar drugs (e.g., two kinase inhibitors with shared scaffolds), verify their embeddings are closer (by cosine similarity or Euclidean distance) than two structurally dissimilar drugs. This is the core functional requirement.
3. **Invalid SMILES handling:** Given an invalid or unparseable SMILES string, verify the module raises a clear error or returns a null/zero embedding (documented behavior, not a silent crash).
4. **Differentiability:** Verify gradients flow through the encoder. Given a loss on drug_emb, verify gradients are computed for the encoder's parameters.
5. **Determinism:** Given the same SMILES, verify the encoder produces the same embedding in eval mode (no stochasticity in forward pass).
6. **Batch independence:** Verify that encoding drug A in a batch with drug B produces the same embedding as encoding drug A alone (no inter-drug interaction in the encoder).

## Integration test

Validate M2 is pretraining-ready and produces structure-aware embeddings:
- Given a panel of ~20 known drugs spanning multiple classes (KIs, anthracyclines, non-cardiac drugs), verify that the embedding space clusters drugs by structural similarity (e.g., KIs with shared scaffolds cluster together).
- Verify that M2 can encode drugs not in the training set (the unseen-drug scenario) — the encoder should produce a valid embedding for any valid SMILES.
- Verify that M2's output dimensionality D is compatible with M3's expected input (M3 takes drug_emb + dose as input).
- Verify that M2 can be frozen (gradient disabled) for Phase 2 fine-tuning scenarios where the encoder is kept fixed.

## Acceptance criteria

- [ ] M2 accepts SMILES strings and produces fixed-dimensional drug_emb
- [ ] Structure-similarity preservation verified (similar drugs → similar embeddings)
- [ ] Invalid SMILES handled gracefully
- [ ] Gradients flow through M2 end-to-end
- [ ] M2 can be frozen for Phase 2
- [ ] Output dimensionality compatible with M3 input

## Risks & notes

- The specific architecture choice (message-passing neural network, transformer, fingerprint-MLP, etc.) is deliberately left open. The key constraint is structure-similarity preservation. The choice should be made based on available code, compute budget, and pretraining data characteristics.
- L1000 pretraining teaches the encoder to map structure → response-relevant features, but L1000 cell lines are cancer cell lines — drug mechanisms may differ in cardiomyocytes. Phase 2 fine-tuning adapts the full model.
- Some cardiomyocyte drugs (e.g., monoclonal antibodies in DToxS-iPSC) don't have small-molecule SMILES. These cannot be encoded by M2 and must be excluded from training or handled with a separate encoding strategy.
- The embedding dimension D is a design decision — large enough to capture chemical diversity, small enough to avoid overfitting with ~2k cardiomyocyte samples.
