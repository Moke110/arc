# AGENTS.md

AI-driven MIRI drug discovery research pipeline: drug SMILES → HG (biological) + UniMol (3D) embeddings → XPert DEG prediction (978 L1000 landmark genes) → MIRI efficacy ranking. Research codebase, not a conventional software project.

## Docs convention

- `AGENTS.md` is the only agent doc, kept at this root (routing index) and in `pert_model/xpert/`, `drug_emb_model/hg/`, `miri_reposition/xpert-fda-v1/` (currently missing — optional cleanup), and `miri_reposition/xpert-fda-v2/`. Read the subdir `AGENTS.md` before touching code there.
- Update any `AGENTS.md` you loaded at session end. No `CLAUDE.md`/`README.md` exist.

## Layout

- `pert_data/dtoxs/` — DToxS dataset (data/docs only, no code). `pert_model/new-model/` — design docs, **not implemented**; don't build on it.
- `drug_emb_model/hg/` — HeteroGNN drug embedding model.
- `pert_model/xpert/` — XPert drug perturbation prediction model (library).
- `miri_reposition/xpert-fda-v1/` — the active MIRI screening workflow (steps 01–11), parity-constrained to the author reference.
- `miri_reposition/xpert-fda-v2/` — the HG-aligned MIRI screening workflow (docs/adr/0006): same 11 steps, but step 02 uses author-pinned frozen inference for the FDA-extension HG rows. `input/` is a symlink to v1's `input/`.
- `disease_data/miri/{GSE}/` — shared per-GSE disease inputs (raw data + DEG lists).
- `miri_reposition/xpert-fda-v1_history/` — frozen pre-refactor pipeline snapshot (gitignored; parity reference).
- `reference/` — literature PDFs, background reading, gitignored.

## Domain glossary

- **author reference** — the L1000 HG embeddings (rows 0–8980) that the author-given XPert was pretrained on: `drug_emb_model/hg/embeddings/article/HG_drug_embeddings.npy`, byte-identical to the article's figshare release. The HG *model* that produced it was never released.
- **FDA extension** — the 2,123 FDA drugs (rows 8981–11103) and their graph inputs (mol-token node features, ChEMBL DTI edges, Tanimoto>0.5 DDS edges).
- **author-pinned inference** — the v2 step-02 frozen layer-by-layer forward over the extended graph, overwriting L1000-drug and gene states with per-layer pins after each convolution so only FDA drug states evolve.
- **HG-aligned workflow** — `miri_reposition/xpert-fda-v2/`; keeps the L1000 side completely frozen and embeds FDA drugs with the same frozen model + same graph recipes.
- **norm mapping** — v1's undocumented scale-patching of FDA rows to the author scale; v2 replaces it with a documented full-CDF quantile alignment of the FDA rows (L1000 rows are never rescaled).
- **parity** — v1's constraint that L1000 rows are byte-identical to the author reference; v2 preserves it by construction (gate B).

## Environment

- Single conda env **`arc`** (Python 3.13) for all subprojects; every `run.sh` calls `conda activate arc`. Created from root `environment.yml` (see `docs/adr/0003`). On this machine it exists as a rename of the legacy `xpert` env (full stack incl. torch+CUDA; R packages for `disease_data/miri/` may be incomplete — install from `environment.yml` if ortholog remapping is needed).
- Subprojects exchange data only via pre-computed files (`.npy`, `.csv`, `.h5ad`); no shared runtime.

## Entrypoints

- `pert_model/xpert/` — `bash pretrain/run.sh` or `python -m src.train`. YAML configs resolve against the owning project root: model-owned pretrain configs relative to `pert_model/xpert/`, workflow-owned finetune configs relative to the workflow.
- `drug_emb_model/hg/` — `bash run.sh` (pretrain); `python -m src.finetune_hg --update-xpert-embeddings --embeddings-output <path>` writes the all-drug HG table (HG never writes into a workflow dir at runtime).
- `miri_reposition/xpert-fda-v1/` — step-wise: 11 step directories (`01_drug_list/`…`11_integrate/`), each with its script, `output/` data, and plots, reading inputs from previous steps + static files in `input/`; orchestrated by `screen.ipynb`. CLI: `python <step>/<step>.py`.
- `miri_reposition/xpert-fda-v2/` — same structure as v1; step 02 (`02_drug_embeddings/02_drug_embeddings.py`) runs the HG-aligned author-pinned inference + `_lib/verify_gates.py` (gates A–D).

## Data facts and gotchas

- Drug embeddings cover 11,104 indices (0–11103). First 8,981 seen in L1000 pretraining; 8,981–11,103 are FDA-extended. In v1 the FDA rows are HG fine-tuned + quantile-norm-mapped; in v2 they are author-pinned inference + full-CDF quantile-aligned (docs/adr/0006).
- The author HG model that produced the author reference was never released (figshare ships only XPert checkpoints + the HG embedding). The saved HG checkpoint (`drug_emb_model/hg/model/20260518_102121/`) reproduces its own phase1 reference, not the author reference.
- Both L1000 and DTOXS are log2 expression space — no z-score conversion in training. The ranking pipeline z-scores the predicted DEG (workflow step 07); always use it.
- The workflow's predict indexes the precomputed UniMol array by drug index — do not regenerate UniMol features from SMILES at predict time (feature drift).
- Known limitation: XPert predicts ~97.2% down-regulated DEGs with very low inter-drug variance, hurting downstream discrimination.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`, `.pdf`) and `.env` are gitignored — never commit them. Pipeline CSVs/PNGs are regenerated outputs.
- No lint/test/typecheck/CI; verification is via notebooks and manual output inspection. Do not invent test commands.
