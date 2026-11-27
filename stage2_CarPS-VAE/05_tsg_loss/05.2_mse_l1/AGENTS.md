# TranSiGen additive Delta L1 sweep

For runner or scoring changes, read the [shared TranSiGen APIs](../../util/README.md#transigen-training-and-evaluation).

Read [README.md](README.md), [Stage 2 routing](../../AGENTS.md), ADRs
0062 and 0086–0090 and the [execution environment](../../../docs/agents/execution-environment.md).

Use the frozen denoise-5k training release, same-plate controls, KPGT vectors
and cold-drug folds from Step 05.1. Select each checkpoint on its registered
validation objective. Run registration, train, predict and replay in separate
processes with one memory-capped CUDA fit at a time. Keep candidates under
the checkout's root `tmp/` and record WSL RAM and GPU usage every ten seconds.

For publication, verify all 27 runs, all three zero-weight controls against
Step 05.1 native MSE, all 243 fold scores and byte-identical CSV-only renders.
