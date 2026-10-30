# AGENTS.md — 03.2_hvg_sensitivity_5k

Routing index. The one-page human summary is `README.md`; the decision is
ADR 0028. This substep is a **sensitivity fork** of 03.1: same wells, folds,
entrants, seeds and budgets — only the gene axis differs.

`execution/`
- `script/runner.py` — fold driver: fit → predict → group-average → evaluate.
- `script/build_gene_set.py` — the 5,000-gene HVG selection + provenance manifest.
- `script/build_data.py` — frozen-input builder; verifies 18k parity before writing.
- `script/compare_18k.py` — paired per-drug/per-fold comparison against 03.1.
- `script/evaluate.py` — AUPRC/PCC/Spearman/R²/RMSE + 2,000-draw bootstrap.
- `script/run.sh` — one model/fold (`--dry-run` = five CPU baselines, fold 0).
- `script/run_queue.sh` — 12 entrants × 5 folds under one GPU flock.
- `script/adapters/<model>/` — one subprocess CLI per entrant; model code is read-only.
- `configs/models.json` — 12-entry registry; `models/` holds author configs + ledgers.
- `tests/` — 53 stdlib tests (adds the gene-set contract to 03.1's suite).
- `Makefile` — `make test`, `make dry-run`, `make gene-set`, `make build-data`, `make compare`.

`output/`
- `input/` — frozen 5k input; `manifest.json` (historical benchmark id
  `temp_env_all_5k`) and `gene_set_5k.manifest.json` are frozen evidence.
- `units/` — per-unit run dirs + shared `_cases/` (regenerable, gitignored).
- `results/` — evaluated tables incl. `comparison_18k.csv`; see `results/README.md`.

`summary.png` / `summary_table.csv` — human figure and table at substep root.

Commands, from the repository root:

```bash
python -m unittest discover -s stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/tests -v
python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/make_summary.py
python stage2/03_models_cmbench/03.2_hvg_sensitivity_5k/execution/script/validate_run.py
```

`repo_root.py` finds the repository root by walking up to `.git`; never add
`Path.parents` indexing or `../..` arithmetic. Absolute AUPRC is **not** comparable with
03.1 — compare paired gaps only.
