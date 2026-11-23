# TranSiGen dose-input optimization — round1

Read [README.md](README.md), [Stage 2 routing](../AGENTS.md) and
[ADR0078](../../docs/adr/0078-dual-source-axis-local-dose-round1.md) before changes.

**Active raw-first comparison.** Use `execution/raw_first.py` for paired preparation,
the tenv KPGT wrapper and fixed-prediction evaluation. The default input is the
published full HVG-DEG view. Render the versioned summary CSV with
`figure/raw_first/render.py`. See the [raw-first workflow](README.md#raw-first-inputs-evaluation-and-figures).
Use the repository resource caps for whole-source preparation. Preserve the accepted comparison records and CSV-only figures referenced in [README.md](README.md).
