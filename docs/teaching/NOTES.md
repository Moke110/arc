# NOTES.md

- Preferred notation for training stages: `[input]+[input]+[model] -> [output]` arrow flowcharts.
- Deliverables: HTML lesson + plain-text ascii reference doc.
- Teaching workspace lives in `docs/teaching/` (subfolder, not repo root).
- User already knows pipeline data facts (11,104 drug indices, 8981 L1000 / 2123 FDA, log2 space, 978 landmarks) — do not re-teach.
- v0 = legacy `workflow/xpert-fda-miri-v0/` + pre-refactor `pert_model/xpert/finetune|predict/` (git `47e1604^`); v1 = finalized step-wise `workflow/xpert-fda-miri-v1/`.
- The pre-refactor snapshot `xpert-fda-v1_history/` is currently absent on disk; v0 scripts + git history are the authoritative v0 source.
