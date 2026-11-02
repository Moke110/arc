# AGENTS.md — 06_prnet_transfer/06.7_cm_tahoe (ticket #7)

Routing index (ADR 0030: ≤30 lines). Purpose, recipe, results and measurements live in
[`README.md`](README.md); the decision is [ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md),
the frozen input is [`06.3_input`](../06.3_input/README.md), the warm-start source is
[`06.5_tahoe_pretrain`](../06.5_tahoe_pretrain/README.md), the sibling random-init matrix is
[`06.6_cm_random`](../06.6_cm_random/README.md), and all training code is
[`06.4_pipeline`](../06.4_pipeline/README.md) (this substep has no trainer of its own).

- `execution/warmstart.py` — the strict warm-start loader: the pipeline's gene-axis/structure
  checks plus `ptlib.checkpoint.load_provenance` (every parameter/buffer key present, shapes
  identical, loaded weights bit-identical), the payload-stage/seed guards and the
  anti-random-init canary; returns the recorded provenance block (sha256, config digest, keys).
- `execution/matrix.py` — pure definition of the matrix (5 folds × 3 seeds × 3 LRs from the
  frozen `locked_config.json`, seed→`#5` checkpoint mapping), unit ids/args, status derivation
  (`ok`/`capped`/`failed`) and the per-fold LR rule (mean of three seeds' best validation
  AUPRC@200; lower LR on a tie).
- `execution/check_init.py` / `build_queue.py` — verify the three #5 checkpoints and write
  `output/units.json` (45 units, each `--init pretrained --init-checkpoint <seed's model.pt>`),
  validated with the pipeline's own `ptlib.queue.load_queue`; fail closed, never a partial load.
- `execution/run_matrix.py` — the queue wrapper: pre-flight, 1 Hz host/GPU peaks, the run-level
  record, the bounded `--profile` mode (staged snapshot + raw axis-check/startup output + the
  #6 path comparison) and `--verify-resume`.
- `execution/write_record.py` — derives `output/record/*` from `output/runs/*` + `units.json`
  and gates it (`--check`): warm-start provenance per candidate/seed, one-to-one config/
  checkpoint/data identity, selection follows from the curves, 15 selected, #6 symmetry,
  no test-path artifacts.
- `execution/tests/` — stdlib `unittest` (67 tests, offline synthetic checkpoints/runs).
- `execution/Makefile`, `execution/run.sh` — `test`, `tables`, `init`, `queue`, `plan`,
  `matrix`, `profile`, `record`, `check`, `resume`, `check-input`.
- `output/` — gitignored compute output; only `profile/` (tracked wrapper evidence) and the
  post-matrix `record/` are force-added.

Rules: never re-implement training/metrics here; never read test treated expression (the test
split belongs to ticket #9); every candidate is reported, no seed is dropped; the formal matrix
must not start before the Lead signals that ticket #5's sweep finished; the heavy run is
`systemd-run … -p MemoryMax=12G` and its peaks go into the record.
