# AGENTS.md — 06_prnet_transfer/eligibility

Routing index for the 460-structure pharmacological eligibility audit (ticket #1):
which temp_env compounds are approved/clinical drugs, mechanism tools, or other
uses. Purpose, the locked rule text, counts and reproduce commands live in
[README.md](README.md); this file only says where things are.

- `execution/build_eligibility.py` — reads the frozen compounds table, collects
  evidence, applies the rule, writes the record in `output/`; `--check` re-runs
  the gates on a written record. `--offline` rebuilds from the local cache.
- `execution/rules.py` — the locked rule (`RULE_TEXT`, `RULE_VERSION`,
  `rule_sha256()`), ATC normalisation and `decide()`. Pure functions, no I/O.
- `execution/evidence.py` — PubChem / Broad Repurposing Hub / BindingDB evidence,
  the response cache (body + JSONL index under `output/cache/`, local),
  identity checks and `frozen_keys()` (467 name records → 460 keys).
- `execution/tests/` — stdlib `unittest`: rule engine, ATC normalisation, identity
  mapping, and the 460-row completeness/consistency gates.
- `execution/run.sh`, `execution/Makefile` — `build`, `offline`, `check`, `test`.
- `output/` — gitignored compute output; the tracked record inside it is
  force-added (`eligibility_table.csv`, `handoff.json`, `summary.json`,
  `old_tier_comparison.csv`, `evidence_manifest.csv`, `manifest.json`).

Rules: every artifact carries the rule sha256 — change the rule only by bumping
`RULE_VERSION` and rebuilding every artifact. `handoff.json` is the interface
ticket #3 consumes (include / uncertain / exclude lists hashed in
`manifest.json`). The frozen ATC/MoA/MeSH strings are annotation under audit, never
eligibility evidence. Do not edit another ticket's subtree.
