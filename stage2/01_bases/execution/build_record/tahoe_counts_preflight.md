# Bases-only Tahoe raw-count release

Latest user scope keeps only `01_bases`. This builder does not restore a deleted
Stage 2 step or historical native corpus. It directly aggregates immutable Tahoe
cell counts onto the frozen 18164 axis, retaining the source Ensembl/token map.
No expression normalization or significance labels occur here.

Public boundary: `prepare_tahoe.py build|verify`, or
`bases.tahoe.build_tahoe` / `verify_tahoe`. Source defaults are read-only Tahoe
expression and metadata; destination is `output/tahoe-counts-v1/`. The CLI checks
the frozen axis digest and uses the shared source-aware mapping resolver. FAM153B
must retain source token61482/ENSG00000289731 at shared index5091.

Counts are summed for stable `(plate,well,cell_line_id)` identities. Shards can
contain two plates. Marker entries, including their upstream -2 expression
sentinel, are dropped before count validation. Every source expression shard is
SHA256-checked against the immutable source manifest; cell counts are conserved.
Production requires 3388 shards,67018 profiles,1400 controls and65618 treated
profiles. The separate training release excludes300 missing-structure profiles
and has65318 treated rows; it is never rewritten by this builder.

Arrow batches are1024cells; one bounded plate accumulator is int64 and is checked
before int32CSR publication. Completed plate parts have atomic payload+source
seals, so interrupted incomplete plates restart without double counting. Stable
lexical cell-line ordering intentionally replaces the historical Python-set row
ordering; downstream replay joins profile/control IDs, not old integer positions.
The final h5ad has only raw integer X, shared-index var with original source
Ensembl/token identities, and well/line/plate/treatment metadata. It is read with
h5py and CSR windows, never anndata backed mode.

Independent `verify` rereads all raw source cells and uses a sparse grouped
reduction rather than production's elementwise accumulation, then compares every
published profile, cell count and raw count vector. It rechecks source and payload
hashes. Full production build/replay acceptance is pending until the running
service and subsequent independent replay complete.

Three public seam tests pass: hand-worked mixed-plate counts and full replay,
sealed-plate resume with changed-source rejection, and ambiguous token rejection.
Real smoke on source shards0/236/1694 verifies7557 profiles,84674 cells and18164
genes through the same build+verify boundary. Final smoke wall37.02seconds,
peakRSS2767120KiB under10GiB cap. These three shards are not full release acceptance.

The separate original62710-axis microbenchmark measured0.52–0.75seconds per
shard, peak~1.1GiB, but does not include full compression/assembly. Source cache
state is not a proven cold-disk measurement. The337.6GB source pass, plate part
compression, final write and independent replay make a1–2hour planning range
plausible, not guaranteed. Real per-plate timings will supersede that estimate.

Run the long build as a named systemd user service with MemoryMax10G,
MemorySwapMax8G, thread caps2 and `/usr/bin/time -v`; keep its fixed code commit,
unit identity, logs and checkpoint progress. Do not modify a live script.
