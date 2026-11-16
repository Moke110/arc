# PRnet native bounded smoke — 2026-11-16

Issue #228; shared input manifest from `log/frozen-inputs/manifest.json`.
No formal matrix results are claimed. Source PRnet f19174bd is imported unchanged;
Sci-Plex5k defaults and the shared Morgan2048 interface are recorded in the API.

TDD: strict-equality/patience boundary tests failed first because the new module
was absent, then passed after shared EarlyStopping implementation. Boundary cases
cover 50th stale epoch, improvement reset, max-epoch cap and nonfinite rejection.
Native output tests use actual 5000-gene PGM, no treated target interface: fixed
seed exact replay preserves training RNG, another seed changes draws, and negative
Normal samples distinguish native inference from the nonnegative decoder mean.
Four initial tests passed, peak RSS918348KiB under MemoryMax10G/SwapMax8G.

Real frozen-input smoke used seed17/fold0, eight rows per partition, two epochs,
CPU, OMP/MKL/OpenBLAS/NumExpr threads2. Each ran alone under coordinator's 10s
RAM/GPU monitor, systemd MemoryMax10G/SwapMax8G and `/usr/bin/time -v`.

| scenario | epochs | best epoch | validation treated MSE | parameter update | fresh-process replay | wall s | peak RSS KiB |
|---|---:|---:|---:|---|---|---:|---:|
| cold-drug | 2 | 2 | 0.9271705746607016 | yes | exact native Normal prediction | 7.73 | 1100752 |
| warm-start | 2 | 2 | 0.9107131868913882 | yes | exact native Normal prediction | 5.39 | 1115988 |

Smoke candidates are preserved outside the tree at `/tmp/arc-prnet-smoke-cold`
and `/tmp/arc-prnet-smoke-warm`; they are disposable bounded evidence, not accepted
model products or matrix entries. Best checkpoint restored before prediction.
Formal runs still require all shared rows, patience50/max1000 and #231 matrix
acceptance. No upstream checkpoint, old mean-output adapter or Stage3 code used.
