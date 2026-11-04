# AGENTS.md — stage2/05_prnet_pretrain

Routing index. Human record: `README.md`. Decisions:
`docs/adr/0031-tahoe-cm-schema-view-and-prnet-pretrain.md` (view contract +
pretrain run) and `docs/adr/0033-stage2-dataset-corpus.md` (this step owns the
view). Trains, and owns the CM-schema view it trains on. Never import
`pert_model/` or `03_models_cmbench/` — both are replicated in-file on purpose.

- `execution/build_pretrain_view.py` — owns the view builder (moved byte-identical
  from the retired `04_tahoe_data`). Reads `02_dataset_corpus/output/corpus/v2/
  tahoe.h5ad`, `.../tables/{genes,compounds}.csv` and `02_dataset_corpus/
  execution/_parts/`; writes `output/pretrain_view/`. Flags: `--check` paths,
  `--qa` gates, `--verify` replays the projection and compares the content hash
  to `476812c0…` (default off; reads parts, writes nothing), `--force` rebuilds
  the gather cache. The view was **moved, not rebuilt** this round.
- `execution/train_prnet.py` — training entry point (`--check` paths/shape/split,
  `--smoke` 2 epochs, `--replot` redraw, default full run); reads `--view`,
  writes `--output`; root by walking to `.git`.
- `execution/_smoke/make_fixture.py` — 300 x 18,164 synthetic view with the
  §1.2–1.4 schema, generated with the `arc` env so `uns['control_mean']` carries
  the production encoding; a test aid, gitignored, not a deliverable.
- `output/` — reruns overwrite; five human artifacts force-added (`config.json`,
  `fit_history.json`, `summary_table.csv`, `summary.png`, `model_card.md`);
  `run.log` + `checkpoint/model.pt` stay local (ticket 05).
- Recipe: same-source copy of the author's PRnet classes + 03.1 adapter loop —
  do not diverge; Gaussian NLL **may be negative**, never clamped.
- Control input is the raw-count group mean, target `log1p_cp10k` — not parity,
  a recorded deviation. Input `output/pretrain_view/tahoe_cm_schema.h5ad`:
  67,018 x 18,164, CSR int32 X + CSR `log1p_cp10k`, obs per spec §1.3.
- Reproduce commands: `README.md`.
