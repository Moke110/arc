# AGENTS.md — 06_prnet_transfer/06.9_evaluation (ticket #9)

Routing index (ADR 0030: ≤30 lines). Purpose, method, results and measurements live in
[`README.md`](README.md); the study design is [`../README.md`](../README.md) and
[`../AGENTS.md`](../AGENTS.md); the decision is
[ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md). The frozen input is
[`06.3_input`](../06.3_input/README.md); the evaluator, predictor and model-card writer are
[`06.4_pipeline`](../06.4_pipeline/README.md); the evaluated checkpoints are
[`06.6_cm_random`](../06.6_cm_random/README.md) and `06.7_cm_tahoe`; the paired baseline
reference is [`06.8_baselines`](../06.8_baselines/README.md).

- `execution/common.py` — roots, checkpoint-root search, the one data-identity recipe
  (same as #6/#7), hashes, atomic JSON.
- `execution/units.py` / `build_units.py` — the 30-unit manifest: reads the two
  `selected_checkpoints.csv` records, re-verifies every checkpoint sha256 and config
  digest, and fails closed unless 5 folds × 3 seeds × 2 inits are complete.
- `execution/score.py` — the driver: `predict.py` then `evaluate.py` per unit (the #4
  CLIs, never a fork), resume-safe, records wall/peak RSS/GPU peak and failures; the
  `test` phase refuses to run without an explicit frozen-selection confirmation.
- `execution/aggregate.py` / `write_record.py` — per-unit, per-fold/seed, per-drug×dose,
  paired-by-drug (bootstrap over drugs), strata and baseline-comparison tables.
- `execution/checks.py` — record gates; `execution/make_cards.py` — 30 model cards;
  `execution/report.py` — renders the README results block from the frozen predictions.
- `execution/tests/` — stdlib `unittest`, offline synthetic records/eval dirs.
- `output/` — gitignored compute output; only `record/`, `model_cards/` and `figures/`
  are force-added.

Rules: never re-implement training, prediction, metrics or model cards; never select a
model, LR, epoch or seed with a test number; every evaluated unit is reported (no seed or
fold dropped); the heavy run is one unit at a time under `systemd-run` and its peaks go
into the record.
