# TranSiGen loss studies

Read [README.md](README.md) and the [Stage 2 routing](../AGENTS.md).

[Step 05.1](05.1_loss_metrics/AGENTS.md) owns the comparison of TranSiGen
Delta loss metrics and coefficients on the fixed denoise-5k axis.

[Step 05.2](05.2_mse_l1/AGENTS.md) owns the additive Delta L1 weight sweep
on the same axis and protocol.

[Step 05.3](05.3_mse-l1_mse-x/AGENTS.md) owns the shared reconstruction MSE
weight sweep at fixed Delta L1 weight 1.75.

[Step 05.4](05.4_mse-l1_pcc/AGENTS.md) and
[Step 05.5](05.5_mse-l1_cos/AGENTS.md) share the bounded Delta similarity
sweep runner and own their PCC and cosine releases respectively.

Keep study evidence under the numbered study's `log/` directory. The two
smoke receipts pin the [shared coordinator capture](05.4_mse-l1_pcc/log/338-coordinator-resources.jsonl)
under Step 05.4 because one monitor covered both studies.

Use the [current figure verifier](README.md#verify-the-current-figures)
for CSV-only replay of the published PNGs. Use the
[pre-refresh revision](README.md#replay-the-original-formal-publications)
to replay the PNG hashes in the original formal receipts.
