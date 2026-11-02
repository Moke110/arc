# AGENTS.md — 06.4_pipeline

Routing index (ADR 0030: ≤30 lines). Purpose, interface and deviations live in
[`README.md`](README.md); the decision is [ADR 0032](../../../docs/adr/0032-prnet-formal-transfer.md).

- `execution/ptlib/` — the **only** training/evaluation code for the study.
  `config`, `data`, `model`, `train`, `predict`, `metrics`, `evaluate`,
  `manifest` (+ `cache`, `tables`, `checkpoint`, `queue`, `modelcard`,
  `runtime`, `h5io`, `paths`).
- `execution/train.py` — one unit: `--stage pretrain|finetune`, `--out`,
  `--fold`, `--seed`, `--init {random,pretrained}` + `--init-checkpoint`,
  `--structure hidden128|latent128|hidden256|<json>`,
  `--train-drugs main|expanded`, `--cache-mode`, `--smoke`, `--check`.
- `execution/predict.py` — target-free deterministic prediction
  (`--stage --role --fold --checkpoint --out`); Δ = prediction − input baseline.
- `execution/evaluate.py` — truth Δ from matched controls; direction AUPRC.
- `execution/run_queue.py` — resume-safe parallel units (`--queue --workers`).
- `execution/make_model_card.py` — regenerate `model_card.md`.
- `execution/smoke_real.py` / `make gpu-smoke` — Phase-B real-data smoke.
- `execution/tools/` — synthetic fixture, frozen-table CSV cache, Tahoe memmap cache.
- `execution/tests/` — stdlib `unittest` (`make test`).

Rules: never import `pert_model/` or the frozen 05/03.1 steps (same-source copies,
do not diverge); never mutate the frozen 06.3 input; Gaussian NLL may be negative
and is never clamped; selection/early stop is validation direction AUPRC@200
(max); prediction never reads treated expression; CMs use c*, not matched controls,
as the model input at validation/test.
