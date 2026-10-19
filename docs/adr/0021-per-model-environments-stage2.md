# Stage 2 models run in per-model conda environments

**Status: accepted**

`environment.yml` declares a single shared runtime (`arc`) for every component and workflow
in the repository. Stage 2 vendors several model families authored in different years
(chemCPA, biolord, CPA, scGen, State, …) whose PyTorch / dependency requirements conflict
with each other and with the `arc` stack (Python 3.13).

**Decision.** Each Stage 2 model runs in its own pinned conda environment named
`arc-stage2-bringup-<model>` (e.g. `arc-stage2-bringup-chemcpa`); the Stage 2 workflow
drives each model through `subprocess`, mirroring how Stage 1 drives `src.train`.

**Why.** Dependency conflicts across model families are near-certain; isolation keeps each
model reproducible and leaves the shared `arc` environment clean for Stage 1's XPert
runtime.

**Consequences.** This is a deliberate exception to the single-runtime convention stated in
`environment.yml`; the exception is scoped to the Stage 2 workflow (`stage2/`). Model source still lives in `pert_model/` (downloaded, `.git`
removed, then untouched); only the environment is per-model.


**Update (2026-10-18).** The old `cm_bench_v1` benchmark is archived at
`stage2/_archive/cm_bench_v1/`; the active bring-up workflow is
`stage2/`. A venv overlay is acceptable where a full conda env would
duplicate a large torch install, as long as the run card records it.

**Update (2026-10-19, ADR 0025).** The active directory is `stage2/01_setup/`
(renamed from `01_bringup/`); biolord and CMonge keep their dedicated venv
overlays, and the other seven models use the `arc-stage2-bringup-<model>`
conda prefix.
