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
