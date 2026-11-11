# AGENTS.md — stage2/01_bases

Routing for the shared Stage 2 foundation (ADR 0061); summary and commands:
[`README.md`](README.md).

- Gene identity: `execution/bases/identity.py`; release build/replay:
  `execution/build_identity.py`; indexed axis: `output/identity-v1/`.
- Before changing source bindings or moving shared assets, read
  `execution/build_record/pre_migration_inventory.md` and ADR 0030.
- Existing setup, native corpora, expression and evaluator entry points remain
  at the compatibility locations listed in README until their migration.
