# Download — XPert

**method: author-provided (manual)** — there is no public download: the XPert
library is author-given and its upstream version is unknown.

## What is here

`source/` holds the upstream-shaped library as vendored:

- `source/src/` — model library: `train.py`, `model_xpert.py`, `model_utils.py`,
  `dataset.py`, `utils.py`, `config.py`, `metrics.py`, `eval_metrics.py`,
  `predict.py`, `vis_gene_attn.py`, `vis_molecule.py`
- `source/pretrain/run.sh` — the pretrain launch recipe
- `source/LICENSE` — MIT, Copyright (c) 2025 Guo Yue

No download command exists; the library was authored/assembled locally and has
no recorded upstream commit (the author log says `version: "unknown"`).

## Local additions

`local/` holds everything added around the library; each entry is declared in
`local_additions[]` of `author_log.json`:

- `AGENTS.md` (routing index), `logs/` (unimol_tools run logs),
  `metadata/` (L1000 gene/drug tables + PPI gene vectors),
  `metrics/` (per-HG evaluation outputs), `reproducing/` (paper figure
  notebooks), `.claude/` and `.vscode/` (local settings).

## Notes

Reference-only for Stage 2 (the preprint X-Pert model is excluded from Stage 2
deployment; this is the Stage 1 XPert implementation the workflows import). All
59 files that were tracked directly in the old `xpert/` tree left the
git index in ticket 04 because they now live under `source/` or `local/`.
