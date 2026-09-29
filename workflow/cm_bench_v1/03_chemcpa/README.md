# 03_chemcpa — chemCPA adapter (ticket 07)

**Status: pending.** Audit-admitted (`docs/entrant_audit.md`): public pretrained
weights (`chemCPA_models.zip`) + CM-finetune recipe. Not yet runnable here — the
checkpoint host (backblaze) is unreachable from this environment and the model
runs in its own pinned env `arc-stage2-v1-chemcpa` (ADR 0021).

Plan: `adapter.py` subclasses `_lib.model_runner.SubprocessAdapter`; `entry.py`
loads the public checkpoint, maps the per-plate baseline transcriptome into
chemCPA's cell-condition embedding, finetunes on the fold's train wells, and
writes `pred.npz`. Smoke prediction on a fixed held-out drug must match the gene
universe shape before the full run.
