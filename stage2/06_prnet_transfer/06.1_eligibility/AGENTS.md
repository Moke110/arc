# AGENTS.md — 06_prnet_transfer/06.1_eligibility

Routing index for the 460-structure pharmacological eligibility audit (ticket #1): which temp_env compounds are approved/clinical drugs, mechanism tools, or other uses. Purpose, the locked rule text, precedence, counts and reproduce commands live in [README.md](README.md); this file only says where things are.

Status: **v1.3** (2026-11-02, ticket #1 review revision). The v1.2 record is frozen in `output/v1.2/`; `output/current.json` is the pointer to the current record.

- `execution/build_eligibility.py` — reads the frozen compounds table, collects evidence, applies the rule, writes the versioned record under `output/<version>/`; `--check` re-runs the gates on the written record; `--offline` rebuilds byte-identically from the local cache.
- `execution/rules.py` — the locked rule (`RULE_TEXT`, `RULE_VERSION`, `PRECEDENCE_DOMAIN`, `rule_sha256()`), ATC normalisation, the human-use marker and other-use helpers, and `decide()` (returns decision, class, precedence label and the per-row evidence basis). Pure functions, no I/O.
- `execution/evidence.py` — PubChem / Broad Repurposing Hub / BindingDB evidence, the response cache (body + JSONL index under `output/cache/`, local), identity checks and `frozen_keys()` (467 name records → 460 keys).
- `execution/tests/` — stdlib `unittest`: rule engine and precedence (32), evidence extraction incl. the bare PubChem Uses statements (16), and the versioned record/gates/change-report checks (18).
- `execution/run.sh`, `execution/Makefile` — `build`, `offline`, `check`, `test`.
- `output/v1.2/` — the locked v1.2 record, byte-identical (rule, table, handoff, summary, old-tier comparison, evidence manifest, manifest).
- `output/v1.3/` — the current record: the same six files plus `change_report.csv` (every row whose decision moved from v1.2, with its evidence).
- `output/current.json` — machine-readable pointer: current artifact version, rule version/sha256, handoff/table/manifest paths and sha256, and the v1.2 predecessor.
- `output/cache/` — local only (gitignored): PubChem bodies and JSONL index, BindingDB match file, ChEMBL probe.

Rules: every artifact carries the rule sha256 — change the rule only by bumping `RULE_VERSION` and rebuilding every artifact into a new `output/<version>/` (never overwrite a frozen record; carry its manifest as the new manifest's `predecessor`). `output/<current>/handoff.json` (schema `prnet-eligibility-handoff-v2`) is the interface ticket #3 consumes; the include / uncertain / exclude lists and the table are hashed in the version's `manifest.json`. The frozen ATC/MoA/MeSH strings are annotation under audit, never eligibility evidence; a mechanism text alone never satisfies eligibility. Do not edit another ticket's subtree.
