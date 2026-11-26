# AGENTS.md — stage1

Stage 1 MIRI drug-repositioning workflows. Each entry is a self-contained
11-step screen or a frozen snapshot; Stage 2 model setup/selection lives in the
top-level `../stage2_CarPS-VAE/` and is not a dependency here.

- `xpert-fda-miri-v1/` — canonical finalized screen (frozen DToxS finetune data).
- `xpert-fda-miri-v1.1/` — DToxS-LINCS testbed sibling.
- `xpert-fda-miri-v2/` — DToxS-primary + DToxS-iPSC-CM finetune sibling.
- `xpert-fda-miri-v3/` — all-CM finetune sibling (current frontier).
- `xpert-fda-miri-v0/` — frozen legacy snapshot (gitignored).

Git history retains the old Stage 2 benchmark at its original `stage2/_archive/cm_bench_v1/` path.
