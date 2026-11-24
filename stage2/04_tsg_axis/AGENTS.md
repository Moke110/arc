# TranSiGen prediction-axis comparison

Read [README.md](README.md), [Stage 2 routing](../AGENTS.md), ADRs 0062
and 0080, and the [execution environment](../../docs/agents/execution-environment.md).
Use the verified Step 03 tenv release, frozen KPGT vectors and published
cold-drug folds. Run one heavy job at a time with the resource caps.

Keep run candidates under the checkout's root `tmp/`. Register all run identities
before fitting. Run train, predict and replay as separate processes. Use the
registered test row identities at the fixed-prediction boundary.

The [single-dose archive](historical/04_model_view/README.md) preserves the
earlier comparison. Keep its sealed files and original provenance available
for historical replay and the shared learner imports in downstream studies.
