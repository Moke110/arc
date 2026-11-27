# TranSiGen prediction-axis comparison

For runner or scoring changes, read the [shared TranSiGen APIs](../util/README.md#transigen-training-and-evaluation).

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
command to replay the scientific report. Build the current comparison from
sealed fits with `execution/build_comparison.py` and render with
`figure/script/render.py` using only the plotted CSVs. Publish the three
canonical PNGs in `figure/` and their plotted values in `figure/data/`. Keep
smoke evaluations outside the formal comparison.

The active TranSiGen adapter and training utilities live in `output/api/tsg_axis/`.
Verify the sealed scores with the pinned evaluator in the README.
