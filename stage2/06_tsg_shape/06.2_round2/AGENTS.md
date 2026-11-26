# TranSiGen prediction architecture comparison

Read [README.md](README.md), the [Stage 2 routing](../../AGENTS.md), ADRs
[0096](../../../docs/adr/0096-transigen-shape-round2.md) and
[0095](../../../docs/adr/0095-transigen-shape-encoder-width-extension.md),
and the [execution environment](../../../docs/agents/execution-environment.md).

Keep registered conditions, matched controls, KPGT and the training settings
fixed. Register each coordinate before training. Train one memory-capped CUDA
coordinate at a time and run prediction and replay in fresh processes.
Record WSL RAM, available host RAM and GPU use at ten-second intervals.
Supply the renderer with verified CSVs.

- Before changing the architecture runtime, reproduce the [bounded CUDA smoke](README.md#reproduce-the-bounded-cuda-smoke).
- Before publishing, verify the [formal matrix and three sealed originals](README.md#reproduce-the-formal-matrix).
- To redraw published figures, use the [plot CSVs and renderer](README.md#replay-the-published-figures).
