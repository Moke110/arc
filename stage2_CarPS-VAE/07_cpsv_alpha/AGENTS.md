# CarPS-VAE alpha

Read the [release and replay instructions](README.md), [ADR 0099](../../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md), and [execution limits](../../docs/agents/execution-environment.md).

Run the public verifier for all nine accepted checkpoints before accepting the release or a study that depends on it. Keep the original Step 06.4 checkpoint directories read-only during replay. A successor registers its chemical representation, fit configuration and evidence in its own step and calls the frozen alpha interface.
