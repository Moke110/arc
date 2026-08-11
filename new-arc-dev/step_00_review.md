# Step 0 — Literature Review & Architecture Decision

## Goal

Systematically review how existing virtual-cell drug- and gene-perturbation models integrate multi-source transcriptome data across platforms, cell types, and batches, and use that evidence to decide the composite architecture for the cardiomyocyte drug-perturbation model.

## Dependencies

None. This is the foundational research step.

## Inputs

- User-provided landscape HTML of 35 drug-perturbation models (2020–2026, A-layer + B-layer)
- Published literature on gene-perturbation foundation models (scGPT, GEARS, scFoundation, CellFlow, Geneformer)
- Three independent benchmarking studies (Ahlmann-Eltze 2025, PerturBench 2024, Wei 2025)

## Outputs

- A catalog of six integration mechanisms with representative models and relevance to this task
- Critical benchmarking evidence on what actually works
- A composite architecture recommendation with evidence basis for each component

## Behavior

### 0.1 Review scope

Review 35 drug-perturbation models from the landscape HTML plus 5 gene-perturbation foundation models and 3 benchmarking studies. For each model, extract: data preprocessing, cell/baseline-state module, drug-encoding module, interaction module, loss function, training strategy, and prior-knowledge integration — focusing on mechanisms that enable cross-platform, cross-cell-type, and cross-batch training.

### 0.2 Six integration mechanisms identified

| Mechanism | What it does | Representative models | Relevance |
|---|---|---|---|
| **A. Adversarial covariate disentanglement** | Encoder decomposes expression into additive basal-state + perturbation + covariate embeddings; adversarial classifiers force basal state to be covariate-invariant | CPA, chemCPA, scCycleMol | Directly addresses platform/batch heterogeneity; CPA is canonical for OOD covariate combos |
| **B. Platform-invariant Δ-embeddings + manifold anchoring** | Maps L1000 and RNA-seq into shared latent space via probabilistic manifold anchoring, preserving platform-specific variance | AetherCell | Directly bridges L1000↔RNA-seq platform gap — our pretrain→finetune scenario |
| **C. Foundation-model pretraining + adapter fine-tuning** | Pretrain transformer on tens of millions of cells, then adapt with drug-conditional adapter | scGPT→scDCA, scFoundation, State, MAP | Transferable cell-state representations; scDCA showed zero-shot to unseen cell lines |
| **D. Cross-attention gene–drug interaction** | Self-attention over genes + cross-attention between drug embeddings and genes | OmniPert, PertDiff, scDCA | Captures cell-type-specific drug mechanisms |
| **E. Knowledge-graph / contrastive pretraining** | Pretrains drug and gene embeddings via contrastive learning over biomedical KG | MAP (187k drugs, 23k genes, 694k relations), SETComp | Improves unseen-drug generalization via mechanistic similarity |
| **F. Flow-matching / diffusion generative modeling** | Models distribution of perturbed states, not just mean | CellFlow, PertDiff, StateXDiff, CaLDM | Autoencoders beat linear on MMD (distributional) metrics |

### 0.3 Critical benchmarking evidence (what actually works)

**Ahlmann-Eltze 2025 (Nature Methods):** Benchmark of GEARS, scGPT, scFoundation, CPA, Geneformer, scBERT, UCE vs. simple linear baselines for gene perturbation prediction.
- No deep-learning model outperformed simple linear baselines on double-perturbation or unseen-single-perturbation tasks.
- Foundation-model pretraining on single-cell atlases provided only a small benefit over random embeddings.
- Pretraining on perturbation data (not atlas data) was the only transfer that consistently helped.
- Implication: pretraining signal must come from perturbation data, not observational atlases.

**PerturBench 2024 (Altos Labs):** Benchmark of CPA, SAMS-VAE, BioLord, GEARS + simple baselines across 6 datasets including covariate-transfer tasks.
- No single model architecture consistently outperforms others. Simpler architectures scale better with larger datasets.
- CPA without its adversarial component (CPA-noAdv) outperformed CPA with it on rank metrics.
- Autoencoder architectures outperform linear baselines on MMD (distributional) metrics but not on average-effect or rank metrics.
- Implication: for cross-cell-line transfer, a simple additive latent model with good covariate embeddings may match complex VAEs; the adversarial component should be ablated, not assumed helpful.

**Wei 2025 (Nature Methods):** Benchmark of 27 methods across 29 datasets with explicit cellular-context-generalization scenarios.
- Cellular context generalization (cross-cell-type transfer) is the hardest scenario; performance correlates negatively with inter-context heterogeneity.
- Methods with explicit cellular context embedding (scDisInFact, biolord, scVIDR) handle multi-condition datasets better.
- Implication: explicit cell-type/context embedding + sufficient target-context fine-tuning data are key for cross-platform cardiomyocyte generalization.

### 0.4 Architecture recommendation

No single published model is optimal. The evidence points to a composite architecture:

**Recommended: chemCPA-style compositional autoencoder + AetherCell-style platform-invariant Δ-embeddings + L1000-pretrained drug encoder + explicit cell-type embedding**

| Component | Source model | Why | Evidence |
|---|---|---|---|
| Compositional latent decomposition (basal state + drug perturbation + covariate, additive) | CPA / chemCPA | Disentangles biology from batch; enables OOD covariate combos; drug encoded by molecular structure for unseen-drug generalization | CPA designed for this; chemCPA extends to unseen drugs |
| Drug molecular encoder (pretrained on L1000 drug space) | chemCPA, MiTCP, MAP | Enables unseen-drug prediction via molecular similarity; L1000 pretrain gives ~20k drug coverage | Ahlmann-Eltze: perturbation-data pretraining helps; MAP-KG contrastive pretraining |
| Platform-invariant Δ-embedding (predict Δ = perturbed − control, not absolute expression) | AetherCell | Bridges L1000 (pretrain) → RNA-seq cardiomyocyte (finetune) platform gap; Δ-space is more platform-invariant | AetherCell; also PRnet predicts fold-changes |
| Cell-type / context embedding (learnable) | scDisInFact, scDCA, State | Explicit context embedding helps cross-cell-type generalization | Wei 2025: this is the key gap |
| Adversarial disentanglement: ABLATE by default | CPA (remove adversary) | PerturBench showed CPA-noAdv > CPA on rank metrics | PerturBench |

**Why not just use one published model:**
- CPA/chemCPA: best-validated compositional design, open code, but no platform-bridging mechanism and adversarial component may hurt.
- AetherCell: best platform-bridging design, but very new (2026), code availability unverified, no public benchmarks yet.
- scGPT/scDCA: foundation-model transfer is attractive but atlas pretraining gives minimal benefit; perturbation-data pretraining is what helps.
- CellFlow: strongest distributional modeling, but designed for genetic perturbations, not drug/chemical.

The composite combines chemCPA's validated compositional structure (with adversary ablated) + AetherCell's Δ-embedding platform bridging + a L1000-pretrained drug encoder + explicit cell-type embedding. This directly targets the two requirements: (1) cross-platform/batch integration via Δ-embeddings + compositional covariate disentanglement, (2) unseen-drug generalization via molecular drug encoding pretrained on L1000's ~20k drugs.

## Unit tests

1. **Model coverage:** Verify all 35 drug-perturbation models from the landscape HTML are classified into one of the six integration mechanisms (or "other"). No model is unclassified.
2. **Evidence grounding:** For each architecture recommendation component, verify there is at least one cited benchmarking finding supporting it. No component is recommended without evidence.
3. **Requirement mapping:** Verify both overarching goals (G1: heterogeneity integration, G2: unseen-drug generalization) are addressed by at least one architecture component each.
4. **Benchmarking completeness:** Verify all three benchmarking studies (Ahlmann-Eltze, PerturBench, Wei) are cited and their key findings extracted.
5. **Adversarial ablation rationale:** Verify the decision to ablate the adversarial component is grounded in PerturBench evidence, not assumed.

## Integration test

Validate that the architecture recommendation is internally consistent:
- Every module in the recommended architecture (M1–M5) maps to at least one integration mechanism from the review.
- The Δ-embedding decision (predict Δ not absolute) is consistent with the platform-bridging requirement.
- The drug encoder pretrain decision is consistent with the Ahlmann-Eltze finding (perturbation-data pretraining, not atlas pretraining).
- The adversarial ablation decision is consistent with the PerturBench finding.
- No recommended component contradicts a benchmarking finding without explicit justification.

## Acceptance criteria

- [ ] Six integration mechanisms cataloged with representative models
- [ ] Three benchmarking studies summarized with key findings and implications
- [ ] Composite architecture recommendation with evidence basis for each component
- [ ] Both goals (G1, G2) explicitly addressed
- [ ] "Why not one model" rationale documented

## Risks & notes

- AetherCell code availability is unverified (2026 preprint). The plan assumes we adapt chemCPA's open codebase and implement the Δ-embedding + platform-bridging components ourselves.
- Given Ahlmann-Eltze and PerturBench results, there is a real risk the model does not beat linear baselines. The plan includes mandatory baseline comparisons (step_12).
- The review focuses on drug-perturbation models; gene-perturbation models are included only for transfer-learning evidence, not as direct architecture candidates.
