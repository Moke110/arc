# AGENTS.md — stage2/06_prnet_transfer/qc (ticket #2)

Routing only. What was decided, what it counts and how to rerun it live in
[`README.md`](README.md); the rule text with its evidence is
[`output/qc_rules.json`](output/qc_rules.json) (locked, hashed).

- `execution/rules.py` — the locked rule set (`RULES`), the two pure decision
  functions, and the reason/flag vocabularies. Thresholds are data here, not code.
- `execution/control_matching.py` — matched-control join keys: temp_env
  `(plate, vehicle=DMSO)`, Tahoe `(plate, cell_line_id)`.
- `execution/build_qc.py` — `--lock` (derive + freeze thresholds), `--build`
  (apply them, write tables), `--verify` (re-hash + count invariants).
- `execution/tests/` — stdlib `unittest`: rule engine, control matching, writer
  determinism, count invariants. `execution/Makefile`, `execution/run.sh`.
- `output/` — gitignored compute output with a `git add -f` tracked record:
  `qc_rules.json`, `qc_summary.json`, `handoff.json`, `manifest.json`, the small
  counts tables. `observations.csv.gz` (71,576 rows) is local-only; its sha256 is
  in `manifest.json`.

Rules: `--lock` must run (and be committed) before `--build`, and `--build`
refuses on any rule-payload drift. No threshold may cite DEG counts, response
strength or model performance — the test suite asserts their absence. Frozen
`02`/`04` outputs and `pert_data/*/source/` are read-only here; the tracked record
is the only thing this subtree writes into git.
