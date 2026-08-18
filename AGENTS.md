# AGENTS.md

AI-driven MIRI drug discovery research pipeline: drug SMILES → HG (biological) + UniMol (3D) embeddings → XPert DEG prediction (978 L1000 landmark genes) → MIRI efficacy ranking. Research codebase, not a conventional software project.

## Docs convention

- `AGENTS.md` is the only agent doc, kept at this root (routing index) and in `pert_model/xpert/`, `drug_emb_model/hg/`, `workflow/xpert-fda-miri-v1/`, and `pert_data/`. Read the subdir `AGENTS.md` before touching code/data there. `pert_data/` is fully gitignored (on-disk only, not tracked) since 2026-08-12 — its `AGENTS.md` lives on disk but is not version-controlled.
- Update any `AGENTS.md` you loaded at session end. No `CLAUDE.md`/`README.md` exist.

## Layout

- `pert_data/` — CM drug-perturbation data lake (see `pert_data/AGENTS.md`): `dtoxs_primary/` (DToxS primary: the released LINCS LDS-1593 well-level package `Data/`/`Metadata/`/…, plus `v1-frozen/` — the byte-identical v1-era package + v1's frozen `dtoxs.h5ad`; the GSE146095/GSE146096 raw GEO deposit deleted 2026-08-18), `dtoxs_ipsc/`, `gse244740/`, `gse262419/`, `gse247541/`, and `l1000-xpert/` (XPert pretrain h5ads). `dtoxs-xpert/` and the `cm/` grouping removed 2026-08-17 (ADR 0011); GSE221595 removed 2026-08-12 (user decision). `pert_model/new-model/` — design docs, **not implemented**; don't build on it.
- `drug_emb_model/hg/` — HeteroGNN drug embedding model (**pure code library**, ADR 0008): `src/` + `run.sh` + `run.ipynb` only. Its author HG reference now lives in the workflow `input/`; the orphaned FDA-extension data/weights are archived at `reference/hg-fda-extension-archive/` (gitignored).
- `pert_model/xpert/` — XPert drug perturbation prediction model (pure library: `src/`, `metadata/`, `metrics/`, `reproducing/`, pretrain `run.sh`). The pretrain weight + YAML live in the workflows (`04_xpert_pretrain/`); `cm` is an additive dataset name (union-CM finetune, treated like `dtoxs`).
- `workflow/xpert-fda-miri-v1/` — the **canonical finalized MIRI screening workflow** (docs/adr/0007): author-aligned, L1000-contained FDA screen (1,004 drugs, `drug_idx < 8981`), 11 steps `01_drug_embed`…`11_report`, workflow-owned `input/` (incl. author HG reference + precomputed UniMol array) and pretrain/finetune YAMLs, authoritative Entrez-verified DToxS→L1000 gene map (`_lib/gene_mapping.py`, `_lib/verify_genes.py`). Frozen DToxS-only finetune data (`input/dtoxs.h5ad`; its fold columns are regenerated in place under the unified fold scheme, ADR 0014).
- `workflow/xpert-fda-miri-v2/` — **union-CM finetune sibling** of v1 (docs/adr/0009): the same author-aligned 1,004-drug screen, but the XPert finetune (steps 05/06) trains on the workflow-built union h5ad (step `03_cm_finetune_data`). Author-aligned since docs/adr/0010: the union follows the author time/dose semantics ({3,6,24} h map, real-concentration dose bins, DToxS family at v1's 24 h slot, 144 h dropped), is per-group quantile-aligned to the L1000 pretrain ([0, 15.002] bins). Since 2026-08-17 the DToxS-primary union arm reads the LINCS LDS-1593 well-level package (docs/adr/0011) — the union is gse262419 + gse244740 + dtoxs_primary (LINCS, 662/45) + dtoxs_ipsc = 7,165 samples. 12 steps `01_drug_embed`…`12_report`; v1-identical except the finetune data + the step-03 composition report (`finetune_report.html`). v1, v1.1, and v2 all benchmark the unified warm-start + cold-drug folds (ADR 0014); v2 is self-contained (no apples-to-apples v1-fold eval).
- `workflow/xpert-fda-miri-v1.1/` — **DToxS-LINCS testbed** (docs/adr/0012): a 12-step mirror of v1 (same 12-step layout as v2) whose only change is the finetune data: step `03_cm_finetune_data` rebuilds the DToxS finetune h5ad (662 samples / 45 drugs) from the LINCS LDS-1593 wells, schema-parity-gated vs v1's frozen `input/dtoxs.h5ad` with an auditable delta report. v1 stays frozen and untouched; v1.1 tests the fresh build + renumbering before any change to v1. Benchmarks the unified warm-start + cold-drug folds (ADR 0014).
- `workflow/xpert-fda-miri-v0/` — frozen legacy pipeline snapshot (gitignored; sole historical reference).
- `disease_data/miri/{GSE}/` — shared per-GSE disease inputs (raw data + DEG lists).
- `reference/` — literature PDFs + the gitignored `hg-fda-extension-archive/`; background reading. `reference/pert_data/` holds the CM dataset papers (PDF + extracted markdown + `publications.md`).
- The removed v1/v2/v3 MIRI workflow variants exist only in git history.

## Domain glossary

- **author reference** — the L1000 HG embeddings (rows 0–8980) that the author-given XPert was pretrained on, workflow-owned at `workflow/xpert-fda-miri-v1/input/HG_drug_embeddings.npy`, byte-identical to the article's figshare release. The HG *model* that produced it was never released.
- **FDA extension** — the 2,123 FDA drugs (rows 8981–11103) and their graph inputs; only used by the removed v1/v2 variants. Preserved for provenance in `reference/hg-fda-extension-archive/`.
- **author-pinned inference** — the removed v2 step-02 frozen layer-by-layer forward over the extended graph (docs/adr/0006, superseded).
- **norm mapping** — the removed v1's undocumented scale-patching of FDA rows to the author scale (superseded by ADR 0007).
- **parity** — the finalized screen's constraint that step-01 HG output is byte-identical to the author reference (`np.array_equal`), and that the refactored workflow reproduces v3's outputs.
- **pretrain unit** — the workflow-owned XPert pretrain weight + `l1000.yaml` (v1: `03_xpert_pretrain/`; v2: `04_xpert_pretrain/`); the L1000 pretrain h5ads live in `pert_data/l1000-xpert/data/`.

## Environment

- Single conda env **`arc`** (Python 3.13) for all subprojects; every `run.sh` calls `conda activate arc`. Created from root `environment.yml` (see `docs/adr/0003`). On this machine it exists as a rename of the legacy `xpert` env (full stack incl. torch+CUDA; R packages for `disease_data/miri/` may be incomplete — install from `environment.yml` if ortholog remapping is needed).
- Subprojects exchange data only via pre-computed files (`.npy`, `.csv`, `.h5ad`); no shared runtime.

## Entrypoints

- `pert_model/xpert/` — `bash pretrain/run.sh` or `python -m src.train`. YAML configs resolve against the owning project root: the workflow-owned pretrain/finetune configs resolve against `workflow/xpert-fda-miri-v1/`.
- `drug_emb_model/hg/` — `bash run.sh` (pretrain; needs the archived `data/` restored from `reference/hg-fda-extension-archive/`). Dormant library.
- `workflow/xpert-fda-miri-v1/` — step-wise: 11 step directories (`01_drug_embed/`…`11_report/`), each with its script(s), `output/` data, and plots, reading inputs from previous steps + static files in `input/`; orchestrated by `screen.ipynb`. CLI: `python <step>/<script>.py`. GPU steps 04/05/06; agentic step 11 (external biomni backend).
- `workflow/xpert-fda-miri-v2/` — same pattern, 12 steps (`01_drug_embed/`…`12_report/`); GPU steps 05/06/07; agentic step 12 (external biomni backend). The finetune/benchmark configs resolve against `workflow/xpert-fda-miri-v2/`.

## Data facts and gotchas

- The finalized MIRI screen covers 1,004 L1000-contained FDA drugs (`drug_idx < 8981`); the FDA-extension rows (8981–11103) are out of scope and never predicted.
- The author HG model that produced the author reference was never released (figshare ships only XPert checkpoints + the HG embedding). The archived HG checkpoint (`reference/hg-fda-extension-archive/model/20260518_102121/`) reproduces its own phase1 reference, not the author reference.
- Both L1000 and the CM datasets are log2 expression space — no z-score conversion in training. The ranking pipeline z-scores the predicted DEG (workflow step 07); always use it.
- The workflow's predict indexes the precomputed UniMol array by drug index — do not regenerate UniMol features from SMILES at predict time (feature drift).
- Known limitation: XPert predicts ~97.2% down-regulated DEGs with very low inter-drug variance, hurting downstream discrimination.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`, `.pdf`), `12_report/cache/` agent JSONs, and `.env` are gitignored — never commit them. Pipeline CSVs/PNGs are regenerated outputs.
- No lint/test/typecheck/CI; verification is via gates and manual output inspection. Do not invent test commands.
