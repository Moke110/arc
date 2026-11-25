# TranSiGen expression shape round one

Read [README.md](README.md), the [Stage 2 routing](../../AGENTS.md),
[ADR 0094](../../../docs/adr/0094-transigen-shape-round1.md) and the
[execution environment](../../../docs/agents/execution-environment.md).
Before editing encoder-width coordinates or reusing prior fits, read
[ADR 0095](../../../docs/adr/0095-transigen-shape-encoder-width-extension.md).

Keep the Step 05.2 weight-1.75 training configuration, including the
objective, input pairing, initialization policy, optimizer, and selection
rules. The native upstream TranSiGen files are read only. The architecture
adapter uses the existing latent numeric dose conditioning in Step 04.

Register before fitting. Run registration, training, prediction, and replay
in separate processes. Use one memory-capped CUDA fit at a time, record WSL
RAM and GPU usage at ten-second intervals, and preserve candidates until
acceptance checks finish. Provide the verified figure CSVs to the renderer.

- Verify or replay the [sealed #347 source](README.md#verify-the-accepted-347-source)
  with its original registrations and receipts.
- Fit and verify the [combined release](README.md#reproduce-the-formal-matrix)
  before replacing its stable paths.
- Regenerate only the [published figures](README.md#replay-the-published-figures)
  with the published CSVs.
