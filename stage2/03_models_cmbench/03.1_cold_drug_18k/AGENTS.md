# AGENTS.md — 03.1_cold_drug_18k

Routing index. The one-page human summary is `README.md`; the contract is
ADR 0027.

`execution/`
- `script/runner.py` — fold driver: fit → predict → group-average → evaluate.
- `script/evaluate.py` — AUPRC/PCC/Spearman/R²/RMSE + 2,000-draw bootstrap.
- `script/build_data.py` — frozen-input builder; `--check` re-verifies hashes.
- `script/run.sh` — one model/fold (`--dry-run` = five CPU baselines, fold 0).
- `script/run_queue.sh` — 12 entrants × 5 folds under one GPU flock.
- `script/adapters/<model>/` — one subprocess CLI per entrant; model code is read-only.
- `configs/models.json` — 12-entry registry; `models/` holds author configs + fidelity ledgers.
- `tests/` — 46 stdlib tests (data contract, metrics, harness, adapters).
- `Makefile` — `make test`, `make dry-run`, `make build-data`, `make fidelity`.

`output/`
- `input/` — frozen benchmark input; do **not** regenerate. `manifest.json` is
  frozen evidence with its historical 2026-10 path strings.
- `units/` — per-unit run dirs + shared `_cases/` (regenerable, gitignored).
- `results/` — evaluated tables; `results/README.md` is the artifact legend.

`summary.png` / `summary_table.csv` — human figure and table at substep root.

Reproduce commands: `README.md` (test, rebuild, validate).

`repo_root.py` finds the repository root by walking up to `.git`; never add
`Path.parents` indexing or `../..` arithmetic. 18k is the primary leaderboard — 03.2 only
changes the gene axis.
