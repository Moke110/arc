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
