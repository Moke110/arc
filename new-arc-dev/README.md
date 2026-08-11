# new_arc — Decoupled Step Plans for Cardiomyocyte Virtual-Cell Drug-Perturbation Model

## Overview

This directory contains **14 decoupled step plans** for building a virtual-cell model that predicts transcriptomic drug-perturbation effects on cardiomyocytes. Each step is independently testable (unit tests + integration test), specifies behavior-level inputs/outputs/dependencies, and can be implemented by reading only that step and its declared dependencies.

The original research plan (systematic review, data inventory, architecture, training plan) is distributed across these step files. The full monolithic plan is preserved at `/mnt/results/execution_trace/PLAN.md` as reference context.

## Two overarching goals

- **G1 — Cross-platform/batch heterogeneity integration:** The model must train jointly on cardiomyocyte data from multiple platforms (3'-DGE, bulk RNA-seq, high-throughput RNA-seq, ScreenSeq), studies, donors, and cell types (primary CM, iPSC-CM).
- **G2 — Unseen-drug generalization on cardiomyocytes:** The model must predict cardiomyocyte transcriptomic responses to drugs it has never seen during training, using only their molecular structure.

## Step list

| Step | File | Title | Dependencies |
|---|---|---|---|
| 0 | `step_00_review.md` | Literature review & architecture decision | none |
| 1 | `step_01_data_inventory.md` | Data inventory & accession validation | 0 |
| 2 | `step_02_retrieve_cm.md` | Data retrieval — cardiomyocyte datasets | 1 |
| 3 | `step_03_retrieve_l1000.md` | Data retrieval — L1000 pretrain corpus | 1 |
| 4 | `step_04_harmonize.md` | Data harmonization & preprocessing | 2, 3 |
| 5 | `step_05_drug_split.md` | Drug overlap analysis & unseen-drug split | 1, 4 |
| 6 | `step_06_impl_m1_m4.md` | Model impl — M1 (cell-state encoder) + M4 (covariate embedder) | 0 |
| 7 | `step_07_impl_m2.md` | Model impl — M2 (drug molecular encoder) | 0 |
| 8 | `step_08_impl_m3.md` | Model impl — M3 (perturbation module) | 7 |
| 9 | `step_09_impl_m5.md` | Model impl — M5 (Δ-decoder) | 0 |
| 10 | `step_10_pretrain_l1000.md` | Phase 1 — L1000 drug-encoder pretrain | 3, 4, 7, 8, 9 |
| 11 | `step_11_train_cm.md` | Phase 2 — cardiomyocyte joint training | 4, 5, 6, 7, 8, 9, 10 |
| 12 | `step_12_evaluation.md` | Evaluation — metrics & baselines | 11 |
| 13 | `step_13_figures_report.md` | Figures & report | 12 |

## Dependency graph

```
Step 0  (Review)
  │
  ▼
Step 1  (Data Inventory)
  ├─────────────────┐
  ▼                 ▼
Step 2 (Retrieve   Step 3 (Retrieve
  CM)                L1000)
  │                 │
  └────────┬────────┘
           ▼
Step 4  (Harmonize)          Step 0 ──► Step 6  (Impl M1+M4)
           │                  Step 0 ──► Step 7  (Impl M2)
           ▼                  Step 0 ──► Step 9  (Impl M5)
Step 5  (Drug Split)                │
           │                  Step 7 ──► Step 8  (Impl M3)
           │                       │
           │     ┌─────────────────┘
           │     │
           ▼     ▼
         Step 10 (L1000 Pretrain)
              │
              ▼
         Step 11 (CM Training) ◄── Step 5, Step 6
              │
              ▼
         Step 12 (Evaluation)
              │
              ▼
         Step 13 (Figures/Report)
```

## Step file convention

Each step file follows this structure:

1. **Goal** — one sentence stating what this step achieves
2. **Dependencies** — which steps must complete first, and why
3. **Inputs** — what this step consumes (behavior-level, no file paths)
4. **Outputs** — what this step produces (behavior-level, no file paths)
5. **Behavior** — what to do, described at behavior level (no function signatures)
6. **Unit tests** — specific test cases with expected behavior (pytest-style descriptions)
7. **Integration test** — end-to-end validation that this step's output meets its schema/properties
8. **Acceptance criteria** — checkable conditions for "done"
9. **Risks & notes** — step-specific risks

**Model module steps (6, 7, 8, 9) use a restricted convention** — they specify only **functional design**: input, output, and objective (what the module must achieve and why). Technical implementation details (specific architectures, layer counts, hidden dimensions, activation functions, dropout rates) are deliberately left open for further decision during implementation. This keeps the plan decoupled from premature architecture commitments while preserving the interface contract between modules.

## Decoupling principles

1. **Each step declares only its direct dependencies.** An engineer can implement step N by reading only step N and its declared dependencies — not the full plan.
2. **Inputs and outputs are behavior-level.** Each step specifies what it consumes and produces at a conceptual level (e.g., "harmonized expression matrix in shared gene space with covariate labels"), not exact file paths or function signatures.
3. **Tests are self-contained per step.** Unit tests validate the step's internal logic; integration tests validate the step's end-to-end output against expected properties. Tests for step N do not depend on steps outside N's dependency list.
4. **Risks are localized.** Each step carries only the risks relevant to it.

## Architecture summary (for reference)

The model is a compositional perturbation autoencoder that predicts the transcriptomic **Δ** (change-from-control) induced by a drug on a cardiomyocyte. Five trainable modules compose as:

```
control_expression ──► [M1: Cell-state encoder] ──► z_basal
drug_SMILES ──► [M2: Drug encoder] ──► drug_emb ──► [M3: Perturbation module] ──► z_pert
dose ──────────────────────────────────────────────────────┘
{cell_type, platform, study, donor} ──► [M4: Covariate embedder] ──► z_cov
                    z = z_basal + z_pert + z_cov  (compositional latent)
                                          [M5: Δ-decoder] ──► predicted_Δ
                        predicted_treated = control_expression + predicted_Δ
```

Key design decisions (evidence basis in step_00):
- Predict Δ, not absolute expression (platform-invariant target)
- MMD regularization instead of adversarial disentanglement (PerturBench: CPA-noAdv > CPA)
- Drug encoder pretrained on L1000 perturbation data (Ahlmann-Eltze: perturbation-data pretraining helps)
- Compositional additive latent (CPA/chemCPA design)
- Explicit cell-type/context embedding (Wei 2025: key gap in current methods)
- Small model size (~7M params) — deliberate for ~2k cardiomyocyte samples
