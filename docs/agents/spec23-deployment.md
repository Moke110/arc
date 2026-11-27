# Spec 23 deployment and WSL cleanup

The materialized checkout was verified on 2026-11-27 in Asia/Shanghai. Local
`main` fast-forwarded from `f0d5bb673435b1a0341301c802ceac5094708a54` through
`dev` to `8a6c43c2cd51e6569f28b39c44b1f83ff3e77d29` for the scientific
verification. The acceptance record was prepared on a separate ticket branch.
The existing unstaged deletion of
`tableA_cardiomyocyte_hypoxia_compound_datasets.csv` remained unstaged throughout.

[Spec #379](https://github.com/Moke110/carps/issues/379) ·
[Ticket #383](https://github.com/Moke110/carps/issues/383) ·
[ADR 0105](../adr/0105-stage-cleanup.md)

## Installed outputs before cleanup

- Stage 1 preserved 286 frozen artifacts, 52 external checkpoint symlink
  targets, both predecessor H5ADs, the frozen phenotype reference and the
  original Annotation seed. The prior v2 tree contained 157 tracked paths and
  322 untracked paths. Materialization excluded eight Python and notebook
  caches and copied 23 additional regular files totaling 74,213,019 bytes.
  All 23 matched their original SHA-256 values. The installed schema gate
  passed on the 958 × 978 union. Rebuilding its 21 official fold columns in
  a scratch copy reproduced the installed H5AD byte for byte.
  [Stage 1 verification](https://github.com/Moke110/carps/blob/9e7c85e38d2ab09c9b509ff8b2f30c67ea07ea85/stage1/verification.md) ·
  [Additional v2 checksums](https://github.com/Moke110/carps/blob/9e7c85e38d2ab09c9b509ff8b2f30c67ea07ea85/stage1/additional_v2_payloads.sha256)
- Stage 2 passed 11 TranSiGen and 25 CarPS-VAE receipt checks, the raw and
  active tenv and Tahoe release checks, and the alpha verifier. Steps 08,
  09 and 10 reproduced 567, 720 and 900 scores respectively. Step 09 also
  verified 33,192 condition rows and Step 10 verified 13,830. The shared
  figure gate reproduced all 45 current PNGs from the sealed tables.
  [TranSiGen verification](../../stage2_CarPS-VAE/util/training/verification.md) ·
  [CarPS-VAE verification](../../stage2_CarPS-VAE/util/training/carps-vae-verification.md)

The Stage 2 input smoke checks for archive absence. The installed pre-cleanup
run temporarily moved the exact history directory within ext4 and restored
the same inode. An earlier installed verification exposed a missing
`EarlyStopping` import in Step 10's formal audit. The corrected revision
passed the entire installed pre-cleanup suite in history isolation before
retirement proceeded.

## Physical cleanup

Following exact-scope confirmations and checks, the five former Stage 1
version directories and `stage2_CarPS-VAE/history` were permanently removed.
The retired v0 reports were unique. The archive inventory identified 4,113
historical paths without a proven replay twin. Explicit confirmation also
covered the independent detached `/tmp/arc-stage2-296-replay` worktree.

Each deletion group was checked for active current-user process references,
mounts and ownership. The active scan examined 222 symlinks across eight
workflow and data roots and found zero targets under retired or temporary
paths. No outbound links in the approved WSL deletion paths and no incoming
links from WSL holdouts were found. Five checkout ticket worktrees were
removed through Git with their branches retained, then all 19 agreed contents
of checkout `tmp/` were cleared with the directory itself retained. The
initial auto-review required separate exact confirmation for the four
additional Stage 1 versions and the 19 checkout scratch entries. These
confirmations were obtained before deletion of those groups.

The fresh WSL inventory after checkout tmp was empty found 155 current-user
ARC/CarPS candidates in `/tmp` and `/var/tmp`. The following 22 paths qualified for removal.

- The explicitly approved detached replay at `212d3e03`, whose tracked tip
  is reachable from `dev`.
- Two clean registered worktrees with retained code. The axis documentation
  worktree had no ignored files. The Tahoe worktree had only generated Python
  bytecode and retained the local branch `impl/stage2-tahoe-284`.
- One documented synthetic verification fixture and 18 generated cache
  directories, comprising 17 Matplotlib font caches and one empty Numba
  cache tree.

The remaining 133 paths include the older `arc-stage2-287`,
`arc-stage2-296` and `arc-tahoe-283` worktrees with unverified ignored
payloads, the `arc-cm-*` research packages and other unresolved historical
scratch. Seventeen stale worktree registrations from the absent earlier
checkout were left in Git. All unresolved paths and artifacts were preserved.

| Measurement on the WSL ext4 filesystem | Available bytes |
| --- | ---: |
| Before copying the retained Stage 1 payloads | 16,245,215,232 |
| After materializing Stage 1 | 10,355,007,488 |
| Immediately before retirement | 10,276,098,048 |
| After removing approved v0 and Stage 2 history | 104,182,841,344 |
| After removing the other approved Stage 1 versions | 117,715,378,176 |
| After clearing checkout tmp | 142,441,578,496 |
| After the 22 WSL temporary removals | 220,529,676,288 |

Available ext4 space increased by 210,253,578,240 bytes from the immediate
pre-retirement measurement, and by 204,284,461,056 bytes net of Stage 1
materialization and the acceptance scratch created before retirement. This
measures free space within WSL. Windows disk and VHDX allocation were outside
the measured scope.

## Detached provenance

Before removal, the two unique baseline proof files were copied byte for
byte to the ignored local evidence directory
`.worktrees/spec23-retained-provenance/`, outside the cleaned roots.
`SHA256SUMS` verified both copies before and after cleanup. The same directory
holds the exact before and remaining WSL inventories and cleanup measurements.

| Preserved file | SHA-256 |
| --- | --- |
| `refs.txt` | `ed212bef27dc56a81b1db97c793273ccfc72e24eff6aa6e501515bd9fb56eef0` |
| `steps01-06-per-path.sha256` | `e6afc96eb6af7cf4d6c1562e5be01bd3a0e9371913f3afb79eb4c8f60dbb7e89` |

## Reproduction after cleanup

The installed verification started with both retired trees and every prior
checkout scratch entry absent. It reran the Stage 1 input, checksum,
checkpoint, schema and official-fold checks. The regenerated union H5AD
again matched exactly. Stage 2 reran the active data smoke, alpha, all three
formal verifiers and shared figures. Their original scientific receipts
matched and all 45 figures passed. A second scan of the 222 retained workflow
and data symlinks again found zero targets under retired or temporary paths.

| Formal study | Reproduced score rows | Scientific receipt SHA-256 |
| --- | ---: | --- |
| Step 08 | 567 | `56b8983678f72adc2deb9a34613561428b3bd77507f44b0ae091b1bc7d40a6f2` |
| Step 09 | 720 | `b8d20c28479c49348a90ef7d6c2071d5c741451fc8f5bd3ee4568de5f09472eb` |
| Step 10 | 900 | `c39a2d6a4b08c6b8c22b2d45fa79e7a72fd28f925301e8b23bec68ebe1dbbd1d` |

The tracked-file check scanned 539 executable files with zero holes.
The data path check found zero unresolved path occurrences and both Stage 2
and Stage 3 path audits passed. The ARC environment's `sed` exited 139
even on `--version`. The shell health gates ran with system GNU sed 4.9,
the available `rg`, and separate ARC Python wrappers for the source
verification dependencies. The first isolated system-only PATH attempt
needed `rg` for the Stage 3 check. Both environment failures were reproduced
and the checks completed without changing repository code or the installed
tool configuration.

The offline source gate verified all 66 external items in 5 minutes
18 seconds. The Stage 1 and Stage 2 scientific checks used the ARC
environment with two numerical threads. Heavy work ran in a 10 GiB
WSL memory scope with an 8 GiB swap cap. The largest process RSS
reported in the post-cleanup scientific checks was 1,655,348 KiB.

The original monitor recorded 521 samples through 18:19:53. Among those
samples, the lowest WSL available memory was 15,652,368 KiB, the lowest
available Windows host memory was 15,101,764 KiB, and the lowest free GPU
memory was 9,170 MiB. One coordinator handoff produced a 351-second sampling
gap. The bounded stdout-only monitor covered cleanup after the file
monitor stopped. Among the 130 collected terminal samples from 18:19:49 to
18:43:02, observed minima were 20,397,400 KiB available in WSL,
14,859,080 KiB available on the Windows host, and 13,495 MiB free on the GPU.

The post-cleanup verification and ticket worktree used fresh checkout
scratch. The raw acceptance logs and timings were retained locally in the
evidence directory before that scratch was cleared.
