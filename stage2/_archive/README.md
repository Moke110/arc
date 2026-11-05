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

`03_models_cmbench/` is the sealed ADR 0027/0028 two-substep comparison, moved
byte-for-byte under ADR 0042. Its old internal path strings and documents remain
historical evidence. To inspect its scripts use this archival prefix; active tenv
work lives at [`../03_modes_tenv/`](../03_modes_tenv/README.md) and imports no
archived code. Preservation hashes are recorded in the active step
`execution/input_record/archive.json`.
