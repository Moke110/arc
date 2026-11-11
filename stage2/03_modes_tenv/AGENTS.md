# AGENTS.md — stage2/03_modes_tenv

Routing for the active tenv comparison (ADR 0042; specification #40).
Read `README.md` for the contract and reproduction commands.

- `execution/script/inputs.py` — build/verify frozen inputs and predictor-safe manifests (#42).
- `execution/script/` — artifact CLI, isolated adapters, metrics and scheduler (#43–#46).
- `execution/tests/` — compact CLI contract tests; formal evidence uses the full gene axis.
- `execution/input_record/` — input build and archival preservation evidence.
- `output/input/` — immutable derived input release; `fold_N/evaluator/` is evaluator-only.
- `discovery/` — candidate evidence and qualification ledger (#41).

Training uses `docs/agents/execution-environment.md` limits and scheduler-owned
telemetry. New code depends on this step and corpus inputs, never archived code.

- Active retrieval scoring and figures use [`../evaluation/README.md`](../evaluation/README.md), with `AUPRC`, `AUPRC_positive`, and `AUPRC_negative` on the full CM significance axis.
