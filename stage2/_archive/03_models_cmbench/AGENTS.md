# AGENTS.md — stage2/03_models_cmbench

Routing index. Purpose, both leaderboards, the full layout and the reproduce
commands live in `README.md`; each substep has its own `AGENTS.md` + `README.md`.

- `03.1_cold_drug_18k/` — the primary 18,164-gene cold-drug benchmark.
- `03.2_hvg_sensitivity_5k/` — the 5,000-gene HVG sensitivity fork of the same 60 units.

Both substeps are self-contained and identically shaped: `execution/{script/,
configs/, models/, tests/, Makefile}` + `output/{input/, units/, results/}` +
`summary.png` + `summary_table.csv`. `script/` holds `runner.py`, `evaluate.py`,
`build_data.py`, `build_gene_set.py` (5k), `compare_18k.py` (5k),
`adapters/<model>/run.py`, the shell drivers, `make_summary.py` and `repo_root.py`.

Rules:

- `repo_root.py` is the only way to find the repository root (walk up to `.git`,
  assert it exists). Never add `Path.parents` indexing or `../..` arithmetic (ADR 0030).
- The harness is duplicated on purpose; fix a bug in **both** substeps.
- Adapters are read-only subprocesses; no benchmark module imports model code.
- `output/input/manifest.json` and `output/results/{run_cards,metrics}.csv` are
  frozen evidence: do not rewrite their strings, do not regenerate them.
- `03.1` is the primary leaderboard; `03.2` changes one axis and must never be
  read as a replacement ranking (ADR 0028).
