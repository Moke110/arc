# TranSiGen reconstruction MSE sweep

Read [README.md](README.md), the [Stage 2 routing](../../AGENTS.md),
ADRs 0062 and 0086–0091, and the [execution environment](../../../docs/agents/execution-environment.md).

Keep the Step 05.2 inputs and runtime fixed. Run registration, training,
prediction and independent replay in separate processes. Select each
checkpoint using the registered validation objective. Keep one
memory-capped CUDA fit active, keep candidates under the checkout's root
`tmp/`, and record WSL RAM and GPU usage every ten seconds.

For formal acceptance require 27 fresh fits, exact replay for each
coordinate and exact agreement of all three unit-weight references with
sealed Step 05.2 runs at Delta L1 weight 1.75. Verify the nine metrics,
eligibility, frozen bandwidths and byte-identical CSV-only rendering.
