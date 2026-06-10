# CLAUDE.md

This file is the **LLM context index** for the `arc` project. It is loaded at the start of every conversation. Its purpose is to route the agent to the correct subdirectory CLAUDE.md based on the task, minimizing unnecessary file reads.

See [README.md](README.md) for the human-readable project documentation.

## CLAUDE.md vs README.md Convention

Every directory in this project maintains two documentation files:

- **CLAUDE.md** — LLM context index. English only. No markdown tables. No diagrams. Concise, machine-parseable structure: Purpose, Structure (file tree), Dependencies, Key Files (file → one-line purpose), Entry Points, Notes. Optimized for minimal token consumption during agent routing.
- **README.md** — Human project documentation. Chinese narrative with English technical terms. Uses tables, ASCII flowcharts, and diagrams for clarity. Explains goals, design rationale, methodology, and background.

**Loading rule**: Root CLAUDE.md is always loaded. Subdirectory CLAUDE.md files are loaded on-demand — the agent reads only the ones relevant to the current task, following index pointers from parent CLAUDE.md files. README.md files are NEVER loaded automatically; they are updated at the end of each conversation session for any directories whose CLAUDE.md was touched.

**Update rule**: At session end, update all CLAUDE.md files that were loaded during the session (to reflect any changes made), then update the corresponding README.md files.

## Project Goal

Predict drugs for myocardial ischemia-reperfusion injury (MIRI) treatment using AI and bioinformatics. Core workflow: predict transcriptional perturbation (978 landmark genes) on cardiomyocytes from small-molecule chemical structures, then rank drugs by MIRI efficacy and cardiotoxicity scores.

## Pipeline

Drug SMILES → UniMol (3D molecular features) + HG (biological embeddings) → XPert (DEG prediction, 978 genes) → fda_miri_scores/ (drug efficacy scoring, 3 pipelines deployed)

## Subdirectories

- [xpert/](xpert/CLAUDE.md) — Drug perturbation prediction model (XPert). Dual-branch Transformer, ~25M params. Status: complete (pretrain + finetune + FDA prediction). Splits: cold-drug (10-fold), cold-cate (6-fold), warm-start (10-fold), dtoxs (all-train). Finetuning: 4 epochs with new HG embeddings.
- [dtoxs/](dtoxs/CLAUDE.md) — Cardiomyocyte drug perturbation dataset (DToxS LDS-1593). 914 samples. Data/docs only, no code.
- [hg/](hg/CLAUDE.md) — Heterogeneous GNN (HeteroGNN). Contrastive learning for drug biological embeddings. Status: pretrained + fine-tuned on extended graph. Full-drug HG embeddings (11,104×256) deployed to xpert/drug_embeddings/.
- [fda_miri_scores/](fda_miri_scores/CLAUDE.md) — MIRI drug efficacy scoring. 3 methods deployed (WCS+RGES, Zhang/ssCMap, Pathway RGES) producing 4 rank components with z-score normalization. Reference gene sets from 3 species. Results: multi-target TKIs (Cabozantinib, Ponatinib) are top candidates. Status: scoring operational.
- [references/](references/CLAUDE.md) — Reference PDF archive. Background reading, no code dependency.

## Cross-Project Dependencies

- hg → xpert/drug_embeddings/all_drugs_idx2HG.npy (biological embeddings, 11104×256)
- UniMol → xpert/drug_embeddings/all_drugs_unimol_arr.npy (3D features, 11104×122×514, generated via unimol_tools)
- xpert → fda_miri_scores/ (3127 drugs × 978 genes DEG predictions)
- GEO → fda_miri_scores/geo/ — GSE160516, GSE108940, GSE137293 raw/processed data
- Features exchanged as pre-computed .npy files — no shared runtime.
- All subprojects use conda env xpert (Python 3.13).
- fda_miri_scores/ R ortholog mapping uses conda env xpert (R + homologene).

## Key Conventions

- Both L1000 and DTOXS data are in log2 expression space — no z-score conversion.
- Drug embeddings cover 11,104 drugs (indices 0–11103). First 8,981 used in L1000 pretraining; 8,981–11,103 are FDA-extended.
- HG embeddings: training rows (0-8980) = author reference; FDA rows = fine-tuned + quantile norm mapping (FDA/Train std=1.000).
- FDA prediction: ensemble of finetuned models, predictions predominantly down-regulated (mean deg ≈ −0.38).

## Current Progress (2026-06-10)

All upstream stages complete; drug efficacy scoring operational:
- HG pretraining + FDA fine-tuning + quantile norm mapping (11,104 drugs, FDA/Train std=1.000, DTI 2.48x)
- XPert L1000 pretraining (8,981 drugs)
- XPert DTOXS finetuning (cold-drug/cold-cate/warm-start/dtoxs splits, 4 epochs, new HG embeddings)
- XPert FDA DEG prediction (3,127 drugs × 978 genes)
- Model evaluation on DTOXS: 4 model groups across 14 deg-only metrics
- GEO MIRI transcriptome data collected (3 species)
- 3 reference gene set JSONs + cross-species consensus + 16 MIRI pathway gene sets
- 3 scoring methods deployed producing 4 rank components with z-score normalization
- Integrated ranking: top drugs are multi-target TKIs (Cabozantinib, Ponatinib, Ripretinib)
- Known limitation: XPert predicts 97.2% down-regulated DEGs with very low inter-drug variance
- Next: validate scoring results against literature/external databases, refine consensus methodology
