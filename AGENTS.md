# AGENTS.md

Routing index for the ARC repo root: AI-driven MIRI drug-repositioning pipeline (SMILES → HG/UniMol embeddings → XPert DEG prediction → efficacy ranking). Research code, not a software project. Details live in each subdirectory's `AGENTS.md`.

## Where to start

1. **Read ADR 0030** — `docs/adr/0030-external-data-provenance-and-step-layout.md`. It is the law for the external-data layout and the Stage 2 step layout.
2. **Read the category `README.md`** for the data/model you need: `pert_data/`, `pert_model/`, `disease_data/`, `drug_emb_model/`. Each category's `AGENTS.md` is the one-line-per-项 routing index.
3. **Read the step `README.md`** for the Stage 2 step you run; every step carries the same members: `AGENTS.md` (routing), `README.md` (human summary), `execution/` (harness + intermediates), `output/` (artifacts later steps read).
4. **Check health** from the repo root (all offline; none needs the network): `script/check_paths.sh`, `script/check_tracked.sh`, `script/verify_all.sh`.

## External data — four provenance categories (27 项)

Every download source is one **项** under ADR 0030: `script/` (downloader + offline verifier), `source/` (untouched download, byte-identical to upstream, gitignored), `local/` (declared local additions, gitignored), `author_log.json` (tracked provenance record), `manifest.sha256` (tracked per-file hashes over `source/`). The machine-readable origin/size/version is each `author_log.json`; the human inventory of every large ignored asset is [`docs/data-inventory.md`](docs/data-inventory.md).

- `pert_data/` — CM drug-perturbation data lake + XPert L1000 pretrain h5ads (8 项: `Tahoe-100M`, `LINCS-LDS-1593`, `DToxS-iPSC`, `GSE244740`, `GSE247541`, `GSE262419`, `XPert-L1000`, `Annotation`); own `AGENTS.md`.
- `pert_model/` — drug-perturbation prediction model library (14 项: nine deployed, five reference-only); own `AGENTS.md`, deployed/reference/excluded split in `README.md`.
- `disease_data/` — shared per-GSE MIRI disease inputs (3 项: `GSE108940` rat, `GSE137293` pig, `GSE160516` mouse); own `AGENTS.md`.
- `drug_emb_model/` — HG HeteroGNN and KPGT drug-embedding libraries (2 项; code libraries, their `source/` trees stay tracked); own `AGENTS.md`.

## Stage 1

`stage1/` — the completed MIRI screens: frozen `xpert-fda-miri-v0/` (legacy snapshot) and the finalized 11-step workflows `xpert-fda-miri-v1/` (canonical), `-v1.1/`, `-v2/`, `-v3/`. Each workflow owns its scripts, prompts, and outputs. Details: `stage1/AGENTS.md`.

## Stage 2

`stage2/` — model setup and selection; stage overview `stage2/README.md`, routing `stage2/AGENTS.md`. Steps (each with `AGENTS.md`, `README.md`, `execution/`, `output/`, human-facing `.png`/`.csv` at step level):

- `01_setup/` — nine deployed non-preprint models load/train on their native assets; run cards under `execution/<model>/`.
- `02_cm_data/` — frozen model-agnostic CM corpus (`output/corpus/v1/`, ADR 0026) and the five `datasets-1` h5ad views.
- `03_models_cmbench/` — the executed 60-unit cold-drug benchmark, split into independently concluded substeps: `03.1_cold_drug_18k/` (primary) and `03.2_hvg_sensitivity_5k/` (HVG sensitivity fork). Each substep carries its own `execution/`, `output/` (`input/`, `units/`, `results/`), `summary.png`, `summary_table.csv`.
- `04_tahoe_data/` — Tahoe-100M pseudobulk pretrain corpus build from the already-downloaded atlas (ADR 0029).
- Non-step directories: `_envs/` (shared per-model benchmark environments), `audit/` (non-computational records), `_archive/` (superseded `cm_bench_v1` benchmark + archived plan).

## Docs, reference, other trees

- `docs/` — decision records ([`docs/adr/README.md`](docs/adr/README.md), 0009–0030 indexed; 0023 is a numbering gap), the AICM roadmap (`roadmap.md`), the ignored-asset inventory (`data-inventory.md`), the wet-lab record (`wetlab-validation-stage1.md`), agent docs, teaching material.
- `GLOSSARY.md` — domain vocabulary, including the external-item / step-layout terms.
- `tableA_cardiomyocyte_hypoxia_compound_datasets.csv` — orphan dataset table, tracked; the same table is in `docs/data-inventory.md`.
- `archived/` — gitignored provenance archive of orphaned HG FDA-extension data/weights.
- `presentation/` — thesis proposal deck (`开题.pptx`) and its background literature review (`my_review.md`).
- `reference/` — background literature library (PDFs, markdowns, index; gitignored).
- `wet_data/` — Stage-1 wet-lab validation raw data (gitignored); record in `docs/wetlab-validation-stage1.md`.
- `logs/` — runtime logs (Tahoe fetch/QA history; ephemeral).
- `environment.yml` — shared conda environment `arc`; `skills-lock.json` — pinned agent-skills lockfile.
- `.agents/`, `.claude/` — agent skills and local config.
- `.scratch/` (issue tracker + specs), `.serena/`, `.vscode/`, `.env`, `.gitignore` — workspace scratch, tooling, and config.

## Agent skills

### Issue tracker

Issues and specs live as local markdown files under `.scratch/<feature-slug>/`. See `docs/agents/issue-tracker.md`.

### Triage labels

Default canonical labels (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`). See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: root `GLOSSARY.md` + `docs/adr/`. See `docs/agents/domain.md`.
