# TranSiGen architecture composition comparison

For runner or scoring changes, read the [shared TranSiGen APIs](../../util/README.md#transigen-training-and-evaluation).

Read [README.md](README.md), the [Stage 2 routing](../../AGENTS.md),
[ADR 0097](../../../docs/adr/0097-transigen-shape-round3.md),
[ADR 0096](../../../docs/adr/0096-transigen-shape-round2.md), and
the [execution environment](../../../docs/agents/execution-environment.md).

Use one memory-capped CUDA fit at a time. Register all coordinates up front
and train, predict, and replay in separate processes. Record available WSL
and host RAM and GPU use at ten-second intervals. Keep the accepted Step 06.1
and 06.2 products read-only. Supply the renderer with verified CSVs.

Verify the bounded CUDA study before launching the formal matrix. Keep the
training seed distinct from the frozen validation and test sampling seeds.
