# AGENTS.md — stage2/07_autoencoder

- [ADR 0052](../../docs/adr/0052-stage2-autoencoder-reconstruction-boundary.md)
  governs the shared axis, split, fitting, and scoring. Step 08 owns perturbation
  prediction and integration.
- [README](README.md) gives the reconstruction CLI and evidence paths. Read the
  [execution envelope](../../docs/agents/execution-environment.md) before real
  corpus access. Use memory-mapped expression and bounded row batches.
- `execution/script/benchmark.py` owns corpus and split verification, fitted
  artifact provenance and the shared scorer. `execution/tests/` exercises the
  workflow on small fixtures.
- Select fitted artifacts on validation only. Reserve real test scoring for
  final publication after `execution/script/compare.py seal` binds all selected
  artifacts to their validation evidence. Commit the seal before test scoring
  and use `compare.py reload` in a separate process for replay.
