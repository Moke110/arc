# AGENTS.md

AI-driven MIRI drug discovery research pipeline: drug SMILES → HG (biological) + UniMol (3D) embeddings → XPert DEG prediction (978 L1000 landmark genes) → MIRI efficacy ranking. Research codebase, not a conventional software project.

## Docs convention

- `AGENTS.md` is the only agent doc, kept at this root (routing index) and in `pert_model/xpert/`, `drug_emb_model/hg/`, and `workflow/xpert-fda-miri-v1/`. Read the subdir `AGENTS.md` before touching code there.
- Update any `AGENTS.md` you loaded at session end. No `CLAUDE.md`/`README.md` exist.

## Layout

- `pert_data/dtoxs/` — DToxS dataset (data/docs only, no code). `pert_data/l1000/` — L1000 XPert pretrain dataset (dataset home). `pert_model/new-model/` — design docs, **not implemented**; don't build on it.
- `drug_emb_model/hg/` — HeteroGNN drug embedding model (**pure code library**, ADR 0008): `src/` + `run.sh` + `run.ipynb` only. Its author HG reference now lives in the workflow `input/`; the orphaned FDA-extension data/weights are archived at `reference/hg-fda-extension-archive/` (gitignored).
- `pert_model/xpert/` — XPert drug perturbation prediction model (pure library: `src/`, `metadata/`, `metrics/`, `reproducing/`, pretrain `run.sh`). The pretrain weight + YAML live in the workflow (`03_xpert_pretrain/`).
- `workflow/xpert-fda-miri-v1/` — the **canonical finalized MIRI screening workflow** (docs/adr/0007): author-aligned, L1000-contained FDA screen (1,004 drugs, `drug_idx < 8981`), 11 steps `01_drug_embed`…`11_report`, workflow-owned `input/` (incl. author HG reference + precomputed UniMol array) and pretrain/finetune YAMLs, authoritative Entrez-verified DToxS→L1000 gene map (`_lib/gene_mapping.py`, `_lib/verify_genes.py`).
- `workflow/xpert-fda-miri-v0/` — frozen legacy pipeline snapshot (gitignored; sole historical reference).
- `disease_data/miri/{GSE}/` — shared per-GSE disease inputs (raw data + DEG lists).
- `reference/` — literature PDFs + the gitignored `hg-fda-extension-archive/`; background reading.
- The removed v1/v2/v3 MIRI workflow variants exist only in git history.

## Domain glossary

- **author reference** — the L1000 HG embeddings (rows 0–8980) that the author-given XPert was pretrained on, workflow-owned at `workflow/xpert-fda-miri-v1/input/HG_drug_embeddings.npy`, byte-identical to the article's figshare release. The HG *model* that produced it was never released.
- **FDA extension** — the 2,123 FDA drugs (rows 8981–11103) and their graph inputs; only used by the removed v1/v2 variants. Preserved for provenance in `reference/hg-fda-extension-archive/`.
- **author-pinned inference** — the removed v2 step-02 frozen layer-by-layer forward over the extended graph (docs/adr/0006, superseded).
- **norm mapping** — the removed v1's undocumented scale-patching of FDA rows to the author scale (superseded by ADR 0007).
- **parity** — the finalized screen's constraint that step-01 HG output is byte-identical to the author reference (`np.array_equal`), and that the refactored workflow reproduces v3's outputs.
- **pretrain unit** — the workflow-owned XPert pretrain weight + `l1000.yaml` (`03_xpert_pretrain/`); the L1000 pretrain h5ads live in `pert_data/l1000/data/`.

## Environment

- Single conda env **`arc`** (Python 3.13) for all subprojects; every `run.sh` calls `conda activate arc`. Created from root `environment.yml` (see `docs/adr/0003`). On this machine it exists as a rename of the legacy `xpert` env (full stack incl. torch+CUDA; R packages for `disease_data/miri/` may be incomplete — install from `environment.yml` if ortholog remapping is needed).
- Subprojects exchange data only via pre-computed files (`.npy`, `.csv`, `.h5ad`); no shared runtime.

## Entrypoints

- `pert_model/xpert/` — `bash pretrain/run.sh` or `python -m src.train`. YAML configs resolve against the owning project root: the workflow-owned pretrain/finetune configs resolve against `workflow/xpert-fda-miri-v1/`.
- `drug_emb_model/hg/` — `bash run.sh` (pretrain; needs the archived `data/` restored from `reference/hg-fda-extension-archive/`). Dormant library.
- `workflow/xpert-fda-miri-v1/` — step-wise: 11 step directories (`01_drug_embed/`…`11_report/`), each with its script(s), `output/` data, and plots, reading inputs from previous steps + static files in `input/`; orchestrated by `screen.ipynb`. CLI: `python <step>/<script>.py`. GPU steps 04/05/06; agentic step 11 (external biomni backend).

## Data facts and gotchas

- The finalized MIRI screen covers 1,004 L1000-contained FDA drugs (`drug_idx < 8981`); the FDA-extension rows (8981–11103) are out of scope and never predicted.
- The author HG model that produced the author reference was never released (figshare ships only XPert checkpoints + the HG embedding). The archived HG checkpoint (`reference/hg-fda-extension-archive/model/20260518_102121/`) reproduces its own phase1 reference, not the author reference.
- Both L1000 and DTOXS are log2 expression space — no z-score conversion in training. The ranking pipeline z-scores the predicted DEG (workflow step 06); always use it.
- The workflow's predict indexes the precomputed UniMol array by drug index — do not regenerate UniMol features from SMILES at predict time (feature drift).
- Known limitation: XPert predicts ~97.2% down-regulated DEGs with very low inter-drug variance, hurting downstream discrimination.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`, `.pdf`), `11_report/cache/` agent JSONs, and `.env` are gitignored — never commit them. Pipeline CSVs/PNGs are regenerated outputs.
- No lint/test/typecheck/CI; verification is via gates and manual output inspection. Do not invent test commands.
