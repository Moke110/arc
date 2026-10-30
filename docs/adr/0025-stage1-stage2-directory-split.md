# Split stage1 and stage2, rename setup, and exclude preprint models

**Status: accepted (2026-10-19)**

**Update (2026-10-30, ADR 0030):** paths in this record were rewritten by the external-data / step-layout restructure; the decision itself is unchanged.

The repository mixed the Stage 1 MIRI screening workflows and the Stage 2 model
selection work under a single `workflow/` tree. Stage 2 also contained a
subdirectory (`workflow/stage2/`) whose depth made every path relative to the
Stage 1 parent. The 2026-10 survey added several preprint models to
`pert_model/`, which made the deployed-model set ambiguous.

**Decision.**

1. `workflow/` becomes `stage1/` and contains only the Stage 1 MIRI screens
   (`xpert-fda-miri-v0`, `-v1`, `-v1.1`, `-v2`, `-v3`).
2. `workflow/stage2/` becomes top-level `stage2/`. The old Stage 2 benchmark
   `workflow/_archive/cm_bench_v1/` moves to `stage2/_archive/cm_bench_v1/`.
3. `stage2/01_bringup/` becomes `stage2/01_setup/`; active Stage 2 documents use
   **setup**, not **bring-up**.
4. Preprint models are excluded from deployment: X-Pert, DORA, BaiZe, and
   PerturbLDM are removed from `pert_model/` and get no setup run card.
5. The setup smoke set is the nine non-preprint, gate-pass models: PRnet, CPA,
   chemCPA, biolord, PerturbNet, scVIDR, Squidiff, CMonge, DECANT. WAVE and ExPO
   are non-preprint but fail the revised gates (no dose / 978-only), so they are
   reference-only.
6. One dedicated environment per model. Existing names
   `arc-stage2-bringup-<model>` are retained to avoid re-installing the four
   working environments; new models use the same prefix. biolord and CMonge keep
   their dedicated `.venv` overlays (allowed by ADR 0021).

**Consequences.** Paths in live docs, run cards, and scripts are updated to
`stage1/` and `stage2/`; historical ADRs keep their original terminology but
their path references are rewritten. `stage2/01_setup` is a smoke stage, not a
benchmark; no model comparison, fold, metric, or CM adapter is defined here.
Preprint literature remains recorded in `pert_model/SURVEY_2026-10-19_v2.md`, but
its code is no longer vendored.

See `stage2/README.md` and `stage2/01_setup/README.md`.
