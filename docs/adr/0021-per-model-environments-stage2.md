# Stage 2 models run in per-model conda environments

**Status: accepted**

`environment.yml` declares a single shared runtime (`arc`) for every component and workflow
in the repository. Stage 2 vendors several model families authored in different years
(chemCPA, biolord, CPA, scGen, State, …) whose PyTorch / dependency requirements conflict
with each other and with the `arc` stack (Python 3.13).

**Decision.** Each Stage 2 model runs in its own pinned conda environment named
`arc-stage2-<version>-<model>` (e.g. `arc-stage2-v1-chemcpa`); the benchmark workflow drives
each model through `subprocess`, mirroring how Stage 1 drives `src.train`. Stage 2 v2 reuses
a v1 environment for models shared between the two versions.

**Why.** Dependency conflicts across model families are near-certain; isolation keeps each
model reproducible and leaves the shared `arc` environment clean for Stage 1's XPert
runtime.

**Consequences.** This is a deliberate exception to the single-runtime convention stated in
`environment.yml`; the exception is scoped to Stage 2 benchmark workflows
(`workflow/cm_bench_*`). Model source still lives in `pert_model/` (downloaded, `.git`
removed, then untouched); only the environment is per-model.
