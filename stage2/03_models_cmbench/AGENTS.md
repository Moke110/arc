# AGENTS.md — stage2/03_models_cmbench

Routing index. Purpose, both leaderboards and the reproduce commands live in
`README.md`; each substep has its own `README.md` + `AGENTS.md`.

```
03_models_cmbench/
├── 03.1_cold_drug_18k/        the primary 18,164-gene benchmark
└── 03.2_hvg_sensitivity_5k/   the 5,000-gene HVG sensitivity fork
```

Each substep is self-contained and has the same shape:

```
<substep>/execution/
├── script/      runner.py (fold driver), evaluate.py (metrics), build_data.py,
│                build_gene_set.py (5k), adapters/<model>/run.py, run.sh,
│                run_queue.sh, gpu_watch.sh, make_summary.py, repo_root.py
├── configs/models.json   12-entry registry (5 baselines + 7 real models)
├── models/<model>/       author config, run card, FIDELITY.md ledger
├── tests/                stdlib unittest contract suite
└── Makefile              `make test`, `make dry-run`, `make build-data`
<substep>/output/
├── input/       frozen benchmark input (manifest.json, arrays, true sets, splits)
├── units/       per-unit run dirs, `_cases/`, run_status.csv (regenerable)
└── results/     leaderboard/metrics/per_drug/bootstrap/run_cards/fidelity
```

Rules:

- `repo_root.py` is the only way to find the repository root (walk up to `.git`,
  assert it exists). Never add `Path.parents` indexing or `../..` arithmetic — see ADR 0030.
- The harness copy is duplicated on purpose; fix a bug in **both** substeps.
- Adapters are read-only subprocesses; no benchmark module imports model code.
- `output/input/manifest.json` and `output/results/{run_cards,metrics}.csv` are
  frozen evidence: do not rewrite their strings, do not regenerate them.
- `03.1` is the primary leaderboard; `03.2` changes one axis and must never be
  read as a replacement ranking (ADR 0028).
