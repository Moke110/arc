# TranSiGen bounded admission smoke

2026-11-16, ticket #229. Inputs: the #227 frozen shared inputs at
`/tmp/arc-tenv5k-inputs`; both scenario fold0, seed17, first eight rows of each
shared train/validation/test partition, two complete epochs, CPU, two OpenMP/MKL
threads. This is not a formal matrix entry or model-selection result.

Both executions used `execution/transigen.py --smoke-rows 8 --max-epochs 2`,
including its fresh Python subprocess replay. Both updated weights, restored
best_epoch2, generated identity-aware well/parent scores and exactly replayed
native sampled treated predictions. `cold.json` and `warm.json` preserve terminal,
configuration, software environment and artifact hashes. Payload candidates remain
at `/tmp/arc-tenv5k-229-smoke-cold` and `/tmp/arc-tenv5k-229-smoke-warm`.

Coordinator-owned single-heavy-task lock, systemd MemoryMax10G/MemorySwapMax8G,
`/usr/bin/time -v`; preflight available RAM41GiB. Cold peak RSS1,373,792KiB,
wall8.34s; warm peak1,373,128KiB, wall5.10s; both exit0. Native prediction test
also passes with poisoned treated-label access, fixed-seed replay, changed-seed
variation and preservation of caller RNG state (peak806,204KiB).

The first test was red with missing model adapter, then green through actual
unmodified upstream model construction/forward. Upstream selective initialization,
VAE loss, latent sampling and singleton training skip are retained. No source,
pretrained weights, dose inputs or Stage3 implementation were changed.
