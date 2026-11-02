# AGENTS.md — 06_prnet_transfer/06.9_evaluation (ticket #9)

Routing index (ADR 0030: ≤30 lines). Purpose, method, results and measurements live in
[`README.md`](README.md); the design is [`../README.md`](../README.md) and
[`../AGENTS.md`](../AGENTS.md); the decision is
[ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md). Inputs: frozen
[`06.3_input`](../06.3_input/README.md); evaluator, predictor and model-card writer
[`06.4_pipeline`](../06.4_pipeline/README.md); evaluated checkpoints
[`06.6_cm_random`](../06.6_cm_random/README.md) / `06.7_cm_tahoe`; paired baseline
[`06.8_baselines`](../06.8_baselines/README.md).

- `execution/common.py` — repo-root walk, recorded-path helpers (`relpath`/`repo_path`),
  checkpoint-root search, the one data-identity recipe (same as #6/#7), hashes, atomic JSON.
- `execution/units.py` / `build_units.py` — the 30-unit manifest: reads the two
  `selected_checkpoints.csv` records, re-verifies every checkpoint sha256 and config
  digest, and fails closed unless 5 folds × 3 seeds × 2 inits are complete.
- `execution/score.py` — `predict.py` then `evaluate.py` per unit (the #4 CLIs, never a
  fork), resume-safe, records wall/peak RSS/GPU peak and failures; `test` refuses to run
  without an explicit frozen-selection confirmation.
- `execution/aggregate.py` / `write_record.py` — per-unit, per-fold/seed, per-drug×dose,
  paired-by-drug (bootstrap over drugs), strata and baseline-comparison tables.
- `execution/checks.py` — record gates; `make_cards.py` — 30 model cards; `report.py` —
  renders the README results block from the frozen predictions.
- `execution/tests/` — stdlib `unittest`, offline synthetic records/eval dirs. `output/` —
  gitignored compute output; only `record/`, `model_cards/` and `figures/` are force-added.

Rules: never re-implement training, prediction, metrics or model cards; never select a
model, LR, epoch or seed with a test number; every evaluated unit is reported (no seed or
fold dropped); every recorded path is repo-relative, never a checkout absolute; the heavy
run is one unit at a time under `systemd-run` and its peaks go into the record.
