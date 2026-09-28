# AGENTS.md

XPert model library for drug-induced transcriptional perturbation prediction. File index only; ownership/provenance lives in ADRs and workflow docs.

- `src/` — model library: `train.py`, `model_xpert.py`, `model_utils.py`, `dataset.py`, `utils.py`, `config.py`, `metrics.py`, `eval_metrics.py`, `predict.py`, `vis_*.py`.
- `pretrain/` — `run.sh`, the pretrain launch recipe.
- `metadata/` — reference data: `l1000_gene_info_978.csv`, `PPI_gene_vector_128d.npy`, `l1000_sdst_drugs_info_8432.csv`.
- `metrics/` — `evaluate_dtoxs.py` plus `new-hg/` and `old-hg/` per-HG metric tables and plots.
- `reproducing/` — figure notebooks: `fig3/`, `fig4/`, `fig5/`, `fig6/`.
- `logs/` — runtime logs (e.g. `unimol_tools_*.log`).
- `LICENSE` — license text.
- `.claude/`, `.vscode/` — local agent/editor settings.
