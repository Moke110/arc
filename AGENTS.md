# AGENTS.md

AI-driven MIRI drug discovery research pipeline: drug SMILES → HG (biological) + UniMol (3D) embeddings → XPert DEG prediction (978 L1000 landmark genes) → MIRI efficacy/cardiotoxicity ranking. This is a research codebase, not a conventional software project.

## Documentation convention (read this first)

- Every directory keeps **two** docs: `CLAUDE.md` (English, concise, agent index) and `README.md` (Chinese, human narrative). Root `CLAUDE.md` is the routing index for all subdirectories.
- Read the subdirectory `CLAUDE.md` before touching code there (`xpert/`, `hg/`, `dtoxs/`, `fda_miri_vs/`). Do not add tables or diagrams to `CLAUDE.md` files.
- At session end, update every `CLAUDE.md` you loaded to reflect changes, then update the matching `README.md`. `README.md` files are never loaded automatically.

## Environment

- All subprojects assume a conda env named **`xpert`** (Python 3.13); every `run.sh` calls `conda activate xpert`. The env is NOT currently installed on this machine (only base/KPGT/biomni_e1/dargus_e1/hermes_test exist) — create it before running scripts.
- `fda_miri_vs/` R ortholog mapping (homologene) also uses conda env `xpert` (R + homologene).
- Subprojects exchange data only via pre-computed files (`.npy`, `.csv`, `.h5ad`); there is no shared runtime between them.

## Entrypoints

- `xpert/` — model. `bash pretrain/run.sh`, `bash finetune/run.sh <pretrained_path>`, `bash predict/predict.sh`, or `python -m src.train`. YAML configs resolve all paths relative to `xpert/`, so run commands from the subdirectory.
- `hg/` — HeteroGNN. `bash run.sh` (pretrain); `python -m src.finetune_hg --update-xpert-embeddings` writes to `xpert/drug_embeddings/`.
- `fda_miri_vs/` — screening is **Jupyter-notebook-driven**: `fda_miri_virtual-screening.ipynb` (6-step pipeline) and `fda_miri_plots.ipynb`. CLI scoring exists: `python scripts/score_{wcs_rges,zhang,pathway_rges}.py --deg-csv xpert_fda-cm_deg.csv --zscore`. gpt-researcher reports (screening steps 4–6) rely on API keys in `.env` + local embedding model under `.agents/skills/gpt-researcher/`.
- `dtoxs/` — dataset only, no code. `new-model/` — design docs for an XPert loss change; **not implemented** and would require full retraining; don't build on it.

## Data facts and gotchas

- Drug embeddings cover 11,104 indices (0–11103). First 8,981 were seen in L1000 pretraining; 8,981–11,103 are FDA-extended (HG fine-tuned + quantile-norm-mapped).
- Both L1000 and DTOXS are log2 expression space — no z-score conversion anywhere in training. But the ranking scripts take `--zscore` on predicted DEG; always use it.
- Known limitation: XPert predicts ~97.2% down-regulated DEGs with very low inter-drug variance, hurting downstream discrimination.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`, `.pdf`) and `.env` are gitignored — never commit them. Pipeline CSVs/PNGs are regenerated outputs.
- No lint/test/typecheck/CI exists; verification is via notebooks and manual output inspection. Do not invent test commands.
