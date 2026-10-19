# Category taxonomy for the repository root

The root is split into four categories — dev-tool files, drug perturbation prediction components (`pert_data/`, `drug_emb_model/`, `pert_model/`), Stage 1 drug repositioning workflows (`stage1/`), and Stage 2 model setup/selection (`stage2/`) — plus `disease_data/` for disease-specific shared inputs. Components and workflows are category containers holding named instances (e.g. `pert_model/xpert/`, `stage1/xpert-fda-miri-v1/`), not renames of the old flat directories. This gives the taxonomy room for future datasets, embedding models, and per-disease workflows without renaming existing ones.

> Amended 2026-10-19 (ADR 0025): `stage2/` is a top-level stage container for model setup and selection, separate from the Stage 1 workflow container `stage1/`.
