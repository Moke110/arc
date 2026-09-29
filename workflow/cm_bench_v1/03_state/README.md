# 03_state — State adapter (ticket 11)

**Status: pending.** Audit-admitted: Arc Institute public pretrained weights
(SE-600M cell model + transition model) with a finetuning path
(`state emb` / `state tx train`). Runs in pinned env `arc-stage2-v1-state`
(ADR 0021).

Plan: `adapter.py` subclasses the runner; `entry.py` embeds baselines with the
SE-600M checkpoint, finetunes the transition model on the fold's train wells, and
predicts held-out drugs. Smoke prediction must match the gene universe shape.
