# AGENTS.md

AI-driven MIRI drug discovery research pipeline: drug SMILES → HG (biological) + UniMol (3D) embeddings → XPert DEG prediction (978 L1000 landmark genes) → MIRI efficacy ranking. Research codebase, not a conventional software project.

## Docs convention

- `AGENTS.md` is the only agent doc, kept at four locations: this root (routing index), `pert_model/xpert/`, `drug_emb_model/hg/`, and `miri_reposition/xpert-fda-v1/`. Read the subdir `AGENTS.md` before touching code there.
- Update any `AGENTS.md` you loaded at session end. No `CLAUDE.md`/`README.md` exist.

## Layout

- `pert_data/dtoxs/` — DToxS dataset (data/docs only, no code). `pert_model/new-model/` — design docs, **not implemented**; don't build on it.
- `drug_emb_model/hg/` — HeteroGNN drug embedding model.
- `pert_model/xpert/` — XPert drug perturbation prediction model (library).
- `miri_reposition/xpert-fda-v1/` — the active MIRI screening workflow (steps 01–11).
- `disease_data/miri/{GSE}/` — shared per-GSE disease inputs (raw data + DEG lists).
- `miri_reposition/xpert-fda-v1_history/` — frozen pre-refactor pipeline snapshot (gitignored; parity reference).
- `reference/` — literature PDFs, background reading, gitignored.

## Environment

- Single conda env **`arc`** (Python 3.13) for all subprojects; every `run.sh` calls `conda activate arc`. Created from root `environment.yml` (see `docs/adr/0003`). On this machine it exists as a rename of the legacy `xpert` env (full stack incl. torch+CUDA; R packages for `disease_data/miri/` may be incomplete — install from `environment.yml` if ortholog remapping is needed).
- Subprojects exchange data only via pre-computed files (`.npy`, `.csv`, `.h5ad`); no shared runtime.

## Entrypoints

- `pert_model/xpert/` — `bash pretrain/run.sh` or `python -m src.train`. YAML configs resolve against the owning project root: model-owned pretrain configs relative to `pert_model/xpert/`, workflow-owned finetune configs relative to the workflow.
- `drug_emb_model/hg/` — `bash run.sh` (pretrain); `python -m src.finetune_hg --update-xpert-embeddings --embeddings-output <path>` writes the all-drug HG table (HG never writes into a workflow dir at runtime).
- `miri_reposition/xpert-fda-v1/` — step-wise: 11 step directories (`01_drug_list/`…`11_integrate/`), each with its script, `output/` data, and plots, reading inputs from previous steps + static files in `input/`; orchestrated by `screen.ipynb`. CLI: `python <step>/<step>.py`.

## Data facts and gotchas

- Drug embeddings cover 11,104 indices (0–11103). First 8,981 seen in L1000 pretraining; 8,981–11,103 are FDA-extended (HG fine-tuned + quantile-norm-mapped).
- Both L1000 and DTOXS are log2 expression space — no z-score conversion in training. The ranking pipeline z-scores the predicted DEG (workflow step 07); always use it.
- The workflow's predict indexes the precomputed UniMol array by drug index — do not regenerate UniMol features from SMILES at predict time (feature drift).
- Known limitation: XPert predicts ~97.2% down-regulated DEGs with very low inter-drug variance, hurting downstream discrimination.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`, `.pdf`) and `.env` are gitignored — never commit them. Pipeline CSVs/PNGs are regenerated outputs.
- No lint/test/typecheck/CI; verification is via notebooks and manual output inspection. Do not invent test commands.
