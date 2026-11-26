# TranSiGen Delta loss metrics

Read [README.md](README.md), the [Stage 2 routing](../../AGENTS.md), ADRs
0062 and 0086–0089, and the [execution environment](../../../docs/agents/execution-environment.md).

Use the verified Step 03 denoise-5k tenv release, frozen KPGT vectors and
cold-drug folds. Keep the Step 04 TranSiGen architecture, matched plate
controls, fixed validation and test seeds, and held-out scoring boundary.
Select each checkpoint on the objective registered for that coordinate.

Register the grid before fitting. Use separate train, predict and replay
processes and one memory-capped CUDA fit at a time. Keep run candidates under
the checkout's root `tmp/`. Use the frozen Step 04 HVG-DEG 5k bandwidths,
metric order, eligibility counts and equal-fold aggregation.

Verify all three native MSE anchors, the complete score ledger and all three
CSV-only renders before replacing the active result set.
