# stage2/_archive — superseded Stage 2 material

Everything here is kept for provenance and is **read-only**; nothing is an active
dependency of Stage 2. `cm_bench_v1/` is the old whole-transcriptome benchmark
(the attempt that preceded the ADR 0027 contract): its 14-model entrant audit was
moved to [`../audit/`](../audit/README.md) and the rest is the frozen record of
what was run and how. `PLAN.md` is the superseded stage plan, archived as-is —
its historical process narrative and layout box describe the stage before the
ADR 0030 restructure. The live version of its still-relevant content (goal,
current selection decision, 978 quarantine, environments, and what each step
produced) lives in [`../README.md`](../README.md). Do not edit an archived file
to match the current tree: it is the record of what was decided at the time.

The later ADR 0027/0028 benchmark was removed on 2026-11-07. Its tracked
records are in Git at `f95fe1e`. Active tenv work lives at
[`../03_modes_tenv/`](../03_modes_tenv/README.md).

[`07_lxpert/`](07_lxpert/README.md) preserves the original ADR 0051 workflow,
checkpoint, figure and treatment smoke records. The archived files retain
their original bytes and recorded paths. To replay the original CLI at those
paths, use Git revision `922db06` with the inputs described in its records.
Active reconstruction lives in [Step 07](../07_autoencoder/README.md).
