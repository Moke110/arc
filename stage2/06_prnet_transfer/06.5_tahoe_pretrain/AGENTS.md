# AGENTS.md — 06_prnet_transfer/06.5_tahoe_pretrain

Routing index (ADR 0030: ≤ 30 lines). Purpose, record schema, measurements and
operator rules live in [`README.md`](README.md); the study decision is
[ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md) and the only
training/evaluation code is [`06.4_pipeline/`](../06.4_pipeline/README.md).

- `execution/run_pretrain.py` — the sweep wrapper: preflight (frozen input,
  derived memmap cache, frozen-table cache, locked recipe, RAM floor, no
  concurrent pretrain), one `06.4 train.py --stage pretrain` per seed
  (sequential, **no recipe override**), then the record + inline check.
- `execution/plan.py` — resolves each unit through the 06.4 `ptlib` config code
  (per-seed config digest; the documented `smoke` / `one-epoch` profile budgets).
  Never re-derives the recipe, never imports `pert_model/`.
- `execution/record.py` — stdlib record writer/reader/checker (schema
  `06.5_tahoe_pretrain/record/1`): input identity, per-unit validation curve,
  resources, three-seed completeness, no-seed-cherry-picking guard.
- `execution/check_record.py` — strict validator for a written record
  (`make check-formal` / `make check-profile`).
- `execution/tests/` — stdlib `unittest` on synthetic fixtures (`make test`).
- `output/profile/` — short-profile record + model card, tracked (`git add -f`).
- `output/formal/` — the three-seed record, curve CSV and model cards, tracked
  once the Lead-run sweep completes; run trees/checkpoints stay local.

Rules: never mutate the frozen 06.3 input or the derived caches; never edit 06.4
code from here; one pretrain at a time (`run.sh` adds the systemd-run memory
cap); no seed may be dropped from the record; the test split never selects.
