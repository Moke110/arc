# TranSiGen prediction-axis comparison

Read [README.md](README.md), [Stage 2 routing](../AGENTS.md), ADRs 0062,
0080 and 0086, and the [execution environment](../../docs/agents/execution-environment.md).
Use the verified Step 03 tenv release, frozen KPGT vectors and published
cold-drug folds. Run one heavy job at a time with the resource caps.

Keep run candidates under the checkout's root `tmp/`. Register all run identities
before fitting. Run train, predict and replay as separate processes. Use the
registered test row identities and each condition's matched plate control at
the fixed-prediction boundary. Select checkpoints by the native validation
objective with 30-epoch patience.

Use `execution/evaluate.py` for released predictions and run its `verify`
command to replay the report. Render with `figure/script/render.py` using the
report summary CSVs. Publish the accepted matched-control figures under
`figure/`. Keep smoke evaluations outside the formal comparison.

Replay Spec 5 scores and figures using the pinned checkout in the README.

The [single-dose archive](historical/04_model_view/README.md) preserves the
earlier comparison. Keep its sealed files and original provenance available
for historical replay and the shared learner imports in downstream studies.
