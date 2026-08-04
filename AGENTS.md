# AGENTS.md

AI-driven MIRI drug discovery research pipeline: drug SMILES → HG (biological) + UniMol (3D) embeddings → XPert DEG prediction (978 L1000 landmark genes) → MIRI efficacy ranking. This is a research codebase, not a conventional software project.

## Documentation convention (read this first)

- Every directory keeps **two** docs: `CLAUDE.md` (English, concise, agent index) and `README.md` (Chinese, human narrative). Root `CLAUDE.md` is the routing index for all subdirectories.
- Read the subdirectory `CLAUDE.md` before touching code there (`pert_model/xpert/`, `drug_emb_model/hg/`, `pert_data/dtoxs/`, `miri_reposition/xpert-fda-v1/`, `disease_data/miri/`). Do not add tables or diagrams to `CLAUDE.md` files.
- At session end, update every `CLAUDE.md` you loaded to reflect changes, then update the matching `README.md`. `README.md` files are never loaded automatically.

## Environment

- All subprojects assume a single conda env named **`arc`** (Python 3.13); every `run.sh` calls `conda activate arc`. It is created from the root `environment.yml` (see `docs/adr/0003-conda-env-arc.md`). The env is NOT currently installed on this machine (only base/KPGT/biomni_e1/dargus_e1/hermes_test exist) — create it before running scripts.
- `disease_data/miri/` R ortholog mapping (homologene) also uses conda env `arc` (R + homologene).
- Subprojects exchange data only via pre-computed files (`.npy`, `.csv`, `.h5ad`); there is no shared runtime between them.

## Entrypoints

- `pert_model/xpert/` — model library. `bash pretrain/run.sh` or `python -m src.train`. YAML configs resolve against the owning project root: model-owned pretrain configs relative to `pert_model/xpert/`, workflow-owned finetune configs relative to the workflow.
- `drug_emb_model/hg/` — HeteroGNN. `bash run.sh` (pretrain); `python -m src.finetune_hg --update-xpert-embeddings --embeddings-output <path>` writes the all-drug HG table (HG never writes into a workflow dir at runtime).
- `miri_reposition/xpert-fda-v1/` — screening is **numbered-stages-driven** (scripts/01-10), orchestrated by `screen.ipynb`. CLI: `python scripts/0{1,6,7,8,9,10}_*.py`. gpt-researcher reports (stage 09) rely on API keys in `.env` + local embedding model under `.agents/skills/gpt-researcher/`.
- `pert_data/dtoxs/` — dataset only, no code. `pert_model/new-model/` — design docs for an XPert loss change; **not implemented** and would require full retraining; don't build on it.
- `disease_data/miri/` — shared per-GSE disease inputs (raw data + DEG lists), consumed by workflow stage 07.
- `miri_reposition/xpert-fda-v1_history/` — frozen pre-refactor pipeline snapshot (gitignored; parity reference).

## Data facts and gotchas

- Drug embeddings cover 11,104 indices (0–11103). First 8,981 were seen in L1000 pretraining; 8,981–11,103 are FDA-extended (HG fine-tuned + quantile-norm-mapped).
- Both L1000 and DTOXS are log2 expression space — no z-score conversion anywhere in training. But the ranking pipeline z-scores the predicted DEG (workflow stage 06); always use it.
- The workflow's predict indexes the precomputed UniMol array by drug index — do not regenerate UniMol features from SMILES at predict time (feature drift).
- Known limitation: XPert predicts ~97.2% down-regulated DEGs with very low inter-drug variance, hurting downstream discrimination.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`, `.pdf`) and `.env` are gitignored — never commit them. Pipeline CSVs/PNGs are regenerated outputs.
- No lint/test/typecheck/CI exists; verification is via notebooks and manual output inspection. Do not invent test commands.
