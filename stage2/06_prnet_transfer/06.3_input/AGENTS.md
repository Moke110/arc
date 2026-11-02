# AGENTS.md — stage2/06_prnet_transfer/06.3_input

Routing index. The frozen contract, counts, reproduce commands and deviations
live in [`README.md`](README.md); the decision is
[ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md) and every constant
lives in `execution/contract.py`. Ticket #3 owns this substep only.

- `execution/build_input.py` — `--check` (sources + handoffs, writes nothing),
  default (build and freeze), `--verify` (offline re-verification),
  `--stub-handoffs` (synthetic handoffs; tests only, never a frozen build).
- `execution/contract.py` — paths, gene axis, handoff schema, locked config and
  hashing. No other file may hard-code a frozen decision.
- `execution/folds.py` / `consensus.py` / `writers.py` — five folds, per-fold
  `c*`, observation assembly and the h5ad writers; pure functions, testable
  without the 16 GB artifacts.
- `execution/tests/` — stdlib `unittest`; `make test` (run in the `arc` env,
  anndata ≥ 0.11; `ARC_DATA_ROOT=/home/chang/arc` only when running inside a
  tracked-files-only worktree).
- `output/` — frozen artifacts: both `.h5ad` are local-only, the small tracked
  record is force-added (`git add -f`).

Hard rules: never write a frozen 02/04/05 output or a `source/` tree; CM keeps
QC-pass rows only while Tahoe keeps every row flagged; `exclude`/`uncertain`
drugs get no fold; a consensus pool draws only on control rows of plates that
carry a **training** drug; no test-split treated expression enters `c*`, the
folds or the training features.
