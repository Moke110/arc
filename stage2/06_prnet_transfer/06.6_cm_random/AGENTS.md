# AGENTS.md — 06_prnet_transfer/06.6_cm_random (ticket #6)

Routing index (ADR 0030: ≤30 lines). Purpose, recipe, results and measurements live in
[`README.md`](README.md); the decision is [ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md),
the frozen input is [`06.3_input`](../06.3_input/README.md), and all training code is
[`06.4_pipeline`](../06.4_pipeline/README.md) (this substep has no trainer of its own).

- `execution/matrix.py` — pure definition of the matrix (5 folds × 3 seeds × 3 LRs derived
  from the frozen `locked_config.json`), unit ids/args, status derivation (`ok`/`capped`/
  `failed`) and the per-fold LR selection rule (mean of three seeds' best validation
  AUPRC@200; lower LR on a tie).
- `execution/build_queue.py` — writes `output/units.json` (45 units) and validates it with
  the pipeline's own `ptlib.queue.load_queue`.
- `execution/run_matrix.py` — the queue wrapper: pre-flight (`MemAvailable ≥ 12 GiB`, no
  second heavy GPU job), invokes `06.4_pipeline/execution/run_queue.py`, samples host/GPU
  peaks at 1 Hz, writes `output/record/matrix_run.json`; `--verify-resume` proves the
  resume semantics.
- `execution/write_record.py` — derives `output/record/*` from `output/runs/*` and gates it
  (`--check`): 45/45 candidates, one-to-one config/checkpoint/data identity, selection
  follows from the recorded curves, 15 selected checkpoints, no test-path artifacts.
- `execution/tests/` — stdlib `unittest` (31 tests, offline synthetic runs and matrices).
- `execution/Makefile`, `execution/run.sh` — `test`, `tables`, `queue`, `plan`, `matrix`,
  `record`, `check`, `resume`, `check-input`.
- `output/` — gitignored compute output; only `record/` (six small files) is force-added.

Rules: never re-implement training/metrics here; never read test treated expression (the
test split belongs to ticket #9); every candidate is reported, no seed or fold is dropped;
the heavy run is `systemd-run … -p MemoryMax=12G` and its peaks go into the record.
