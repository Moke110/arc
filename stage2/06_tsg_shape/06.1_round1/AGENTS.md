# TranSiGen expression shape round one

Read [README.md](README.md), the [Stage 2 routing](../../AGENTS.md),
[ADR 0094](../../../docs/adr/0094-transigen-shape-round1.md), and the
[execution environment](../../../docs/agents/execution-environment.md).

Keep the Step 05.2 weight-1.75 training configuration, including the
objective, input pairing, initialization policy, optimizer, and selection
rules. The native upstream TranSiGen files are read only. The architecture
adapter uses the existing latent numeric dose conditioning in Step 04.

Register before fitting. Run registration, training, prediction, and replay
in separate processes. Use one memory-capped CUDA fit at a time, record WSL
RAM and GPU usage at ten-second intervals, and preserve candidates until
acceptance checks finish. Provide the verified figure CSVs to the renderer.

The formal driver in `execution/matrix_shape.py` verifies the complete
preregistration before fitting. Use the verified receipt and stable local
`output/formal/` paths for follow-up reads. Recompute scores and replay the
two figures with `execution/report_shape.py verify` against the sealed
Step 05.2 reference. Publish only the scope described in ADR 0094.
