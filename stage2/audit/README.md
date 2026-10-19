# Stage 2 entrant audit

The 14-model selection record, moved from the archived `cm_bench_v1` workflow.
It is **not** a benchmark runner and not a Stage 2 evaluation.

- `entrant_audit.py` — the audit table, capability gates, and Markdown/CSV writer.
- `audit_entrants.py` — regenerate `entrant_audit.md` and `entrant_audit.csv`:
  `python stage2/audit/audit_entrants.py`
- `entrant_audit.md` / `.csv` — the current generated record.

The chemCPA row was corrected after setup: every public checkpoint is an
L1000 **977- or 2,000-gene** `ComPert` state dict (188 drugs, 3 covariates), not
a 5,000/20,570-gene model. chemCPA remains admitted as a Tier S model, but its
whole-transcriptome use requires retraining, not just loading the public zip.
