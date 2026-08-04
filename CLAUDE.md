# CLAUDE.md

This file is the **LLM context index** for the `arc` project. It is loaded at the start of every conversation. Its purpose is to route the agent to the correct subdirectory CLAUDE.md based on the task, minimizing unnecessary file reads.

See [README.md](README.md) for the human-readable project documentation.

## CLAUDE.md vs README.md Convention

Every directory in this project maintains two documentation files:

- **CLAUDE.md** — LLM context index. English only. No markdown tables. No diagrams. Concise, machine-parseable structure: Purpose, Structure (file tree), Dependencies, Key Files (file → one-line purpose), Entry Points, Notes. Optimized for minimal token consumption during agent routing.
- **README.md** — Human project documentation. Chinese narrative with English technical terms. Uses tables, ASCII flowcharts, and diagrams for clarity. Explains goals, design rationale, methodology, and background.

**Loading rule**: Root CLAUDE.md is always loaded. Subdirectory CLAUDE.md files are loaded on-demand — the agent reads only the ones relevant to the current task, following index pointers from parent CLAUDE.md files. README.md files are NEVER loaded automatically; they are updated at the end of each conversation session for any directories whose CLAUDE.md was touched.

**Update rule**: At session end, update all CLAUDE.md files that were loaded during the session (to reflect any changes made), then update the corresponding README.md files.

## Agent skills

### Issue tracker

Issues and specs are tracked as local markdown under `.scratch/<feature>/`. See `docs/agents/issue-tracker.md`.

### Triage labels

Five canonical roles, each label string equal to its name: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context — one `CONTEXT.md` + `docs/adr/` at the repo root. See `docs/agents/domain.md`.

## Project Goal

Predict drugs for myocardial ischemia-reperfusion injury (MIRI) treatment using AI and bioinformatics. Core workflow: predict transcriptional perturbation (978 landmark genes) on cardiomyocytes from small-molecule chemical structures, then rank drugs by MIRI efficacy and cardiotoxicity scores.

## Pipeline

Drug SMILES → UniMol (3D molecular features) + HG (biological embeddings) → XPert (DEG prediction, 978 genes) → fda_miri_vs/ (drug efficacy scoring, 3 pipelines deployed)

## Subdirectories

- [xpert/](xpert/CLAUDE.md) — Drug perturbation prediction model (XPert). Dual-branch Transformer, ~25M params. Status: complete (pretrain + finetune + FDA prediction). Splits: cold-drug (10-fold), cold-cate (6-fold), warm-start (10-fold), dtoxs (all-train). Finetuning: 4 epochs with new HG embeddings.
- [dtoxs/](dtoxs/CLAUDE.md) — Cardiomyocyte drug perturbation dataset (DToxS LDS-1593). 914 samples. Data/docs only, no code.
- [hg/](hg/CLAUDE.md) — Heterogeneous GNN (HeteroGNN). Contrastive learning for drug biological embeddings. Status: pretrained + fine-tuned on extended graph. Full-drug HG embeddings (11,104×256) deployed to xpert/drug_embeddings/.
- [fda_miri_vs/](fda_miri_vs/CLAUDE.md) — MIRI drug screening pipeline. 3 scoring methods (WCS+RGES, Zhang/ssCMap, Pathway RGES) → 4 Z-score components → filtering (WCS<0, critical pathways) → integrated ranking. Results in `fda_screened.csv`. Status: screening operational.
- [references/](references/CLAUDE.md) — Reference PDF archive. Background reading, no code dependency.

## Cross-Project Dependencies

- hg → xpert/drug_embeddings/all_drugs_idx2HG.npy (biological embeddings, 11104×256)
- UniMol → xpert/drug_embeddings/all_drugs_unimol_arr.npy (3D features, 11104×122×514, generated via unimol_tools)
- xpert → fda_miri_vs/ (3127 drugs × 978 genes DEG predictions)
- GEO → fda_miri_vs/geo/ — GSE160516, GSE108940, GSE137293 raw/processed data
- Features exchanged as pre-computed .npy files — no shared runtime.
- All subprojects use conda env xpert (Python 3.13).
- fda_miri_vs/ R ortholog mapping uses conda env xpert (R + homologene).

## Key Conventions

- Both L1000 and DTOXS data are in log2 expression space — no z-score conversion.
- Drug embeddings cover 11,104 drugs (indices 0–11103). First 8,981 used in L1000 pretraining; 8,981–11,103 are FDA-extended.
- HG embeddings: training rows (0-8980) = author reference; FDA rows = fine-tuned + quantile norm mapping (FDA/Train std=1.000).
- FDA prediction: ensemble of finetuned models, predictions predominantly down-regulated (mean deg ≈ −0.38).

## Current Progress (2026-06-16)

All upstream stages complete; drug screening pipeline operational:
- HG pretraining + FDA fine-tuning + quantile norm mapping (11,104 drugs, FDA/Train std=1.000, DTI 2.48x)
- XPert L1000 pretraining (8,981 drugs)
- XPert DTOXS finetuning (cold-drug/cold-cate/warm-start/dtoxs splits, 4 epochs, new HG embeddings)
- XPert FDA DEG prediction (3,127 drugs × 978 genes)
- Model evaluation on DTOXS: 4 model groups across 14 deg-only metrics
- GEO MIRI transcriptome data collected (3 species)
- 3 reference gene set JSONs + cross-species consensus + 19 MIRI pathway gene sets (11 UP + 8 DOWN, 8 critical)
- 3 scoring methods deployed producing 4 Z-score components (WCS_zs, RGES_zs, Zhang_zs, pathway_reversal_zs)
- Screening pipeline (fda_screen.ipynb): filter (WCS<0 + ≥7/8 critical pathways) → 70/3127 drugs pass → rank (weighted integrated Z-score)
- Top candidates: Terazosin, Zolpidem Tartrate, Cefiderocol Sulfate Tosylate, Loracarbef, Erythromycin Glucepate
- Known limitation: XPert predicts 97.2% down-regulated DEGs with very low inter-drug variance
- Next: validate screening results against literature/external databases
