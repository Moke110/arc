# TranSiGen Delta loss comparison

Read [README.md](README.md), [Stage 2 routing](../AGENTS.md), ADRs 0062,
0086–0088 and the [execution environment](../../docs/agents/execution-environment.md).

Use the verified Step 03 denoise-5k tenv release, frozen KPGT vectors and
cold-drug folds. Keep the Step 04 TranSiGen architecture, matched plate
controls, fixed validation and test seeds, and held-out scoring boundary.
Select each checkpoint on the objective registered for that coordinate.

Register the grid before fitting. Use separate train, predict and replay
processes and one memory-capped CUDA fit at a time. Keep run candidates under
the checkout's root `tmp/`. Exclude the bounded smoke fits from the formal
comparison. Score against the frozen Step 04 HVG-DEG 5k bandwidths, preserving
the metric order, eligibility counts and equal-fold aggregation. Render reports
from the accepted CSVs.
