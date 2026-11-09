# AGENTS.md — stage2/07_lxpert

- **Input release** — use `execution/script/lxpert.py prepare` and `verify`
  for the Step 02 raw-count projection, target masks and plate memberships.
  The commands and output contract live in [README.md](README.md).
- **Data access** — read [the execution envelope](../../docs/agents/execution-environment.md)
  before opening real corpora. Use HDF5 row slices for raw counts and memory-map
  the normalized release.
- **Design** — read [ADR 0051](../../docs/adr/0051-latent-xpert.md) for the
  reconstruction and perturbation smoke boundaries.
- **Input tests** — `execution/tests/test_inputs.py` exercises the CLI with
  bounded dense and CSR fixtures.
- **Autoencoder training** — use `fit` and `reload` for the bounded probe.
  Run `freeze`, `develop`, `refit` and `reload` as separate processes for
  formal training and release. The [README](README.md) gives commands and
  artifact paths. Run formal fits serially through the scheduler.
- **Treatment smoke** — run `smoke` against the released refit checkpoint and
  verified Step 03 inputs. Schedule the real CUDA probe. The bounded CPU test
  is `execution/tests/test_treatment_smoke.py`.
