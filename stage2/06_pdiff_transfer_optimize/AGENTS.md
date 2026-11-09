# AGENTS.md — Step 06

- **Replay** the frozen cold-drug comparison through [the baseline commands](README.md#validation-baseline).
- **Interventions and attribution** follow [ADR 0048](../../docs/adr/0048-perturbdiff-transfer-diagnosis.md) and [spec #89](https://github.com/Moke110/arc/issues/89).
- **Schedule** real-data verification, inference and fitting through the dedicated scheduler under the [WSL resource limits](../../docs/agents/execution-environment.md).
- **Publish** result PNGs directly in `output/` and new checkpoints with numbered, readable names in `output/checkpoints/`. Reference reused checkpoints by path and digest.
