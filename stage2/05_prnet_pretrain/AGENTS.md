# AGENTS.md — stage2/05_prnet_pretrain

Routing index. The human record is `README.md`; the contract is
`.scratch/stage2-05-prnet-pretrain/spec.md` §2. This step **trains only**; it
consumes 04's view and rebuilds no data. Never import `pert_model/` or
`03_models_cmbench/` — both are replicated in-file on purpose.

- `execution/train_prnet.py` — the only entry point. `--check` = paths/shape/split
  then exit; `--smoke` = 2 epochs; default = full run. Reads only `--view`,
  writes only `--output`; the repo root is found by walking up to `.git`.
- `execution/run.sh` — one-command reproduction (`PY=` overrides the interpreter).
- `execution/_smoke/make_fixture.py` — a 300 x 18,164 synthetic view with the
  §1.2–1.4 schema, for proving the pipeline before 04's view exists. Generate it
  with the `arc` env (anndata 0.12) so `uns['control_mean']` carries 04's
  `dataframe/0.2.0` encoding; the pinned 0.10.8 reader handles it. Test aid, not
  a deliverable (`*.h5ad` is gitignored).
- `output/` — `run.log`, `fit_history.json`, `summary.png`, `summary_table.csv`,
  `config.json`, `model_card.md`, `checkpoint/model.pt`; reruns overwrite. The
  human-facing five (`config.json`, `fit_history.json`, `summary_table.csv`,
  `summary.png`, `model_card.md`) are force-added to version control; `run.log` and
  `checkpoint/model.pt` stay local (ticket 05).
- Input: `stage2/04_tahoe_data/output/pretrain_view/tahoe_cm_schema.h5ad` —
  67,018 x 18,164, CSR `X` int32 + CSR `layers['log1p_cp10k']`, obs per §1.3,
  `uns['control_mean']` per §1.4; batch = (control mean, rFCFP4 x log10 dose,
  log1p_cp10k row).
- Recipe is a **same-source copy** of `pert_model/PRnet/source/{models/PRnet.py,
  trainer/PRnetTrainer.py,data/_utils.py}` and the 03.1 adapter loop — do not
  diverge; the recipe stays identical, but the reported losses are not commensurable
  with 03.1's. Loss is Gaussian NLL and **may be negative**;
  never clamp, `abs()` or "fix" it.

Reproduce commands: `README.md` (`--check`, `--smoke`, full run).
