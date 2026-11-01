# AGENTS.md — stage2/05_prnet_pretrain

Routing index. The human record is `README.md`; the accepted decision record is
`docs/adr/0031-tahoe-cm-schema-view-and-prnet-pretrain.md` (executed 2026-11-01:
68 epochs, best valid Gaussian NLL −2.628080; contract planned in
`.scratch/stage2-05-prnet-pretrain/spec.md` §2). Trains only; consumes 04's view,
rebuilds no data. Never import `pert_model/` or `03_models_cmbench/` — both are
replicated in-file on purpose.

- `execution/train_prnet.py` — the only entry point. `--check` = paths/shape/split
  then exit; `--smoke` = 2 epochs; `--replot` = redraw from `fit_history.json`;
  default = full run. Reads `--view`, writes `--output`; root by walking to `.git`.
- `execution/run.sh` — one-command reproduction (`PY=` overrides the interpreter).
- `execution/_smoke/make_fixture.py` — a 300 x 18,164 synthetic view with the
  §1.2–1.4 schema, generated with the `arc` env so `uns['control_mean']` carries
  04's `dataframe/0.2.0` encoding; the pinned 0.10.8 reader handles it. Test aid,
  not a deliverable (`*.h5ad` is gitignored).
- `output/` — reruns overwrite; five human-facing artifacts force-added
  (`config.json`, `fit_history.json`, `summary_table.csv`, `summary.png`,
  `model_card.md`); `run.log` and `checkpoint/model.pt` stay local (ticket 05).
- Input: `stage2/04_tahoe_data/output/pretrain_view/tahoe_cm_schema.h5ad` —
  67,018 x 18,164, CSR `X` int32 + CSR `layers['log1p_cp10k']`, obs per §1.3,
  `uns['control_mean']` per §1.4.
- Recipe is a **same-source copy** of `pert_model/PRnet/source/{models/PRnet.py,
  trainer/PRnetTrainer.py,data/_utils.py}` and the 03.1 adapter loop — do not
  diverge. Loss is Gaussian NLL, **may be negative**, never clamped/"fixed".
- Control input is the raw-count group mean, the target `log1p_cp10k` — **not** the
  author's same-matrix parity; a recorded deviation, resolve before any finetune.

Reproduce commands: `README.md` (`--check`, `--smoke`, `--replot`, full run).
