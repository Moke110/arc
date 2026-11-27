# Step 10 CarPS-VAE dose-input comparison

Read the [study design](README.md), [Stage 2 routing](../AGENTS.md),
[ADR 0103](../../docs/adr/0103-carps-vae-dose-inputs.md) and the
[execution limits](../../docs/agents/execution-environment.md).
Use the frozen Denoise-DEG 5k cold-dose-drug folds, matched-plate controls, and
the ten registered dose arms. Keep the Step 07 model release and any sealed
comparison evidence unchanged. Run the bounded CUDA gate before the formal
matrix. Keep the coordinator sampler live from preflight until after the
last fit measurement. Verify the report against the frozen data and replay
figures from their CSVs before accepting a publication.

For model, fit or replay changes, use the [shared APIs](../util/README.md#transigen-training-and-evaluation).
