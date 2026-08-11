# Step 8 — Model Impl: M3 (Perturbation Module)

## Goal

Implement the perturbation module (M3) that transforms the drug molecular embedding + dose into a perturbation vector in the same latent space as the basal state — the bridge between molecular structure and transcriptomic effect.

## Dependencies

- **Step 7 (Impl M2):** M3 consumes M2's output (drug_emb) as its primary input. The interface contract (drug_emb dimensionality D) must be established in step_07 before M3 can be implemented.

## Inputs

- `drug_emb`: D-dimensional drug molecular embedding from M2
- `dose`: scalar (log-transformed concentration)

## Outputs

- `z_pert`: L-dimensional perturbation vector, additively combined with z_basal and z_cov

## Behavior

### M3 — Perturbation module

**Functional objective:** Transform the drug molecular embedding + dose into a perturbation vector in the same latent space as the basal state, representing the transcriptomic change the drug induces.

**Input → Output:** [drug_emb (D-dim); dose (scalar)] → z_pert (L-dim)

**Key requirements:**
- Must combine drug_emb and dose into a single perturbation vector. The dose conditioning allows the model to separate dose effects from drug identity effects — important when integrating datasets with different dose regimes.
- Must output a vector in the same latent space as z_basal (L-dimensional) so they can be additively composed: z = z_basal + z_pert + z_cov.
- No standalone loss — trained end-to-end through the decoder's reconstruction loss. The gradient from Δ-reconstruction flows back through z_pert into M3 and M2.
- The specific architecture (MLP, attention, etc.) is left open for further decision. The key constraint is that it must map drug_emb + dose → a perturbation vector that, when added to a basal state, produces the correct Δ via the decoder.

**Pretraining (Phase 1):** M3 is pretrained jointly with M2 on L1000. During Phase 1, M3 learns to map drug_emb + dose → perturbation vector that reconstructs L1000 Δ-profiles across ~70 cell lines. This teaches it the general mapping from molecular structure + dose → transcriptomic change. During Phase 2, M3 is fine-tuned to adapt this mapping to cardiomyocyte-specific drug responses.

**Contribution to goals:**
- G2 (unseen-drug): M3 is the bridge between molecular structure (M2) and transcriptomic effect. If M2 produces a meaningful embedding for an unseen drug, M3 must map it to a perturbation vector that produces the correct Δ. The L1000 pretraining gives M3 broad exposure to drug→perturbation mappings; cardiomyocyte fine-tuning adapts the mapping to cardiac biology.
- G1 (heterogeneity): The perturbation vector is dose-conditioned, allowing the model to separate dose effects from drug identity — important when integrating datasets with different dose regimes (DToxS clinical concentrations vs. GSE221595 multi-dose vs. Rosell-Hidalgo toxicity-screening concentrations).

## Unit tests

1. **Forward shape:** Given drug_emb (D-dim) and dose (scalar), verify z_pert output has shape (L-dim). Given a batch of B drugs, verify output has shape (B, L).
2. **Dose conditioning:** Given the same drug_emb with two different doses, verify z_pert differs. Verify that dose=0 produces a perturbation vector closer to zero than a non-zero dose (no perturbation at zero dose).
3. **Drug sensitivity:** Given two different drug_embs with the same dose, verify z_pert differs. The module must be sensitive to drug identity, not just dose.
4. **Latent space compatibility:** Verify z_pert has the same dimensionality (L) as z_basal (from M1) and z_cov (from M4), so additive composition is valid.
5. **Differentiability:** Verify gradients flow from z_pert back to drug_emb (into M2) and to M3's own parameters.
6. **Batch independence:** Verify that processing drug A with dose X in a batch with drug B produces the same z_pert as processing drug A alone.

## Integration test

Validate M3 works with M2 in the full composition:
- Given a drug SMILES → M2 → drug_emb → M3 (with dose) → z_pert, verify the full pipeline produces a valid perturbation vector.
- Verify z_pert + z_basal + z_cov has the correct shape for decoder (M5) input.
- Verify that for a held-out drug (not in training), the M2→M3 pipeline still produces a valid z_pert (the unseen-drug scenario).
- Verify that dose scaling works: increasing dose for the same drug changes z_pert in a smooth (non-jumpy) manner.

## Acceptance criteria

- [ ] M3 accepts drug_emb + dose and produces z_pert (L-dim)
- [ ] Dose conditioning verified (dose=0 → minimal perturbation)
- [ ] Drug sensitivity verified (different drugs → different z_pert)
- [ ] z_pert dimensionality matches z_basal and z_cov for additive composition
- [ ] Gradients flow through M3 to M2
- [ ] Full M2→M3 pipeline works for held-out drugs

## Risks & notes

- The dose conditioning design is left open — it could be concatenation, multiplicative gating, or a more complex interaction. The key constraint is that dose=0 should produce minimal perturbation.
- Some cardiomyocyte datasets have a single dose per drug (e.g., DToxS 48h at one concentration), while others have multi-dose (GSE221595). The module must handle both — for single-dose datasets, dose is a constant and the module still benefits from dose conditioning on multi-dose L1000 pretraining.
- The perturbation vector's magnitude should scale with drug effect strength — strongly perturbing drugs should produce larger z_pert vectors than weakly perturbing drugs. This is an emergent property from training, not a hard constraint, but worth monitoring.
