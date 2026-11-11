# #157 gates and verification record

Ticket #157, parent #154. Branch `codex/stage2-bases-157` based on `dev`
`7496743688a9872f8588ad122216b671c4bb45a8`. Worktree
`/home/chang/arc/.scratch/worktrees/stage2-bases-157`; the 34 GB payload tree
stayed in the main checkout `/home/chang/arc` and every data-reading command took
explicit paths. Nothing was copied into the worktree, pushed or merged.

Environment: Python 3.13.13 (`/home/chang/miniconda3/bin/python3`), matplotlib
3.10.9, pandas 2.3.3, NumPy 2.4.6; R 4.5.3 / DESeq2 1.50.2 only as a sealed
reference (no refit ran). Threads pinned to 2 for every command. The documented
`systemd-run --user --scope` form is still unreachable in this session
(`Failed to connect to bus: No data available`), so `execution/capped.sh` applied
the equivalent 10 GiB `RLIMIT_AS` locally for `verify_all.sh`; the measured peak
(77.5 MiB) never approached it.

## A. Final acceptance audit

| id | check | outcome | evidence |
|---|---|---|---|
| A1 | every retained formal product is accepted | **PASS** | the regenerated metadata audit's 44 binding checks resolve each product manifest, envelope, acceptance record and seal digest pair; no accepted product lost a binding when the three payloads were removed |
| A2 | no half-built release remains | **PASS** | `output/` holds exactly the ten stable entries (nine products + `output/api`); the superseded half-built truth candidate was removed (`157-removals.json` R2); `find` finds no `*.building`, `*-v1`, `*-v2` or `build_record` directory anywhere under the step |
| A3 | every figure reproduces independently from its CSV | **PASS** | `figure/data/diagnostic_tenv_pr.csv` (`ff3a82e8…`) → `figure/diagnostic_tenv_pr.png` byte-identical `35d14266…`; [`health/157_figure_reproduction.log`](health/157_figure_reproduction.log) |
| A4 | no unresolved live reference requires a deleted version path | **PASS after finding F1** | the old-path scan over the branch tree (`git grep -nE` on `identity-v1\|native-v1\|shared-v[12]\|training-v[12]\|truth-tenv-v[12]\|tahoe-counts-v1\|hvg-tenv-v1\|hvg-metrics-diagnostic-v1\|build_record`, excluding `log/records`, `log/health`, `log/historical`, `log/155-*`, `log/156-*` and the one-shot migration script) now yields only the audit's own `OLD_NAMES` check list, the *retained* historical mapping directory name, and descriptive prose. F1 fixed the last real defect |
| A5 | required scientific dependencies survive or have verified replacements | **PASS** | the accepted truth-fit dependency on `output/historical-training-release/tenv.h5ad` (`a95bba66…`) and its manifest (`3543df76…`) are byte-identical after the removal; the accepted successor envelope binds only live paths (11/11 bindings resolve); both training corpora and the shared projection are untouched |
| A6 | Stage 1 and the immutable external source trees are unchanged | **PASS** | `git diff --name-only a53c9ff codex/stage2-bases-157 -- pert_data stage1` → 0 files; `script/check_paths.sh` PASS with 0 unresolved literals |
| A7 | health checks completed with resource limits and expression-scan scope respected | **PASS** | section C |

### The regenerated acceptance record (A8)

`log/156-metadata-binding-audit.json` had been captured from the #156 fix
worktree and carried five stale fields (absolute worktree roots and the two
`step_root: true` duplicate flags that only the split checkout produced). It was
regenerated from the main checkout with default single-checkout roots:

```bash
python3 .scratch/worktrees/stage2-bases-157/stage2/01_bases/execution/audit_metadata_bindings.py \
    --step-root /home/chang/arc/stage2/01_bases \
    --write .scratch/stage2-bases/runtime/v3/audit-run<N>.json
```

- `step_root` = `/home/chang/arc/stage2/01_bases`; **no** `.scratch/worktrees`
  string survives anywhere in the record; every absolute `/home/...` path is a
  canonical main-checkout path or the pinned `Rscript`.
- **44 binding checks + 21 structural checks, `non_pass: []`**, and exactly one
  documented `EXPECTED_GAP`: `cross: historical-training-release.inputs.shared_manifest
  -> deleted shared-v1 seal`. It stays a documented historical gap — the payload
  and seal were deleted with the historical shared release before this ticket,
  only the mapping CSVs survive, and the check is not fabricated into a pass.
- Two consecutive runs are **byte-identical**: both
  `sha256 1b5642cba9077fdf04d1173596140bc56d9e5cb4c43fc49c845a72e3389b5344`
  (`cmp` clean), 0.33 s each, peak RSS 113,080 / 113,100 KiB, 0 expression
  payloads opened.
- Seven fields differ from the committed #156 record: the five stale ones plus
  the two that the removal required (the rejected-Tahoe check's name and
  verdict). The check count is unchanged.

## B. Removals

Three payloads totalling **15,274,607,787 B (15.27 GB)** were removed after the
proof required by [`157-removals.json`](157-removals.json):

| id | path | bytes | digest |
|---|---|---:|---|
| R1 | `output/historical-training-release/tahoe_bulk.h5ad` | 11,201,616,617 | `6e3ac2ec…` (pre-removal sha256 = the manifest's declared digest) |
| R2 | `execution/checkpoints/truth-tenv-building` (567 files) | 1,206,080,986 | tree `34d7afc0…` |
| R3 | `execution/checkpoints/tahoe-raw-counts` (28 files) | 2,825,100,184 | tree `715dd384…` |

Disk: 869 G → 855 G used, 275 G → 289 G available. Each removal records the old
path, size, digest, the superseding/rejection evidence with its own sha256, the
consumer scan and the retained evidence. No empty old-path directory needed
removing: the migration had already `rmdir`'d `output/native-v1` and every other
old name was a rename, not a copy.

**Kept deliberately** (recorded in `157-removals.json → evaluated_and_kept`):
`output/shared-expression/` (gate D4 is not satisfied, so retirement is not
authorized), `log/diagnostic/cohorts.csv` (accepted audit evidence; #156
interpretation 4 allowed dropping it, but #157's scope is rejected/superseded
payloads, so it is reported, not deleted) and the now-empty
`execution/checkpoints/` boundary directory.

## C. Health checks

Run one at a time; nothing else heavy was active (`free -g` preflight: 42 GiB
available).

| id | command | scope | outcome | measured |
|---|---|---|---|---|
| C1 | `script/check_paths.sh` (**once**, `GIT_INDEX_FILE` filled by `git read-tree codex/stage2-bases-157`) | branch tree's 939 tracked files; literals under `pert_data/`, `pert_model/`, `disease_data/`, `drug_emb_model/` resolved against `/home/chang/arc` | **PASS**, 0 unresolved (0 distinct); 30 whitelisted occurrences / 20 distinct; 55 stage2 executable files scanned, 0 stale stage2 literals | 8.47 s wall, 6,472 KiB peak RSS — [`health/157_check_paths.log`](health/157_check_paths.log) |
| C2 | `script/check_tracked.sh` (worktree) | every executable file in the branch worktree | **PASS** — 242 code files on disk, 284 tracked, 0 untracked, 0 holes | 0.01 s, 4,604 KiB — [`health/157_check_tracked.log`](health/157_check_tracked.log) |
| C3 | `execution/capped.sh /usr/bin/time -v script/verify_all.sh` (**once**, alone) | the **34 immutable external items** across `pert_data pert_model disease_data drug_emb_model pheno_data`: offline `manifest.sha256` verification of the external items only | **PASS — verified 34/34** | 4:48.32 wall (174.29 s user + 108.15 s sys), 79,376 KiB peak RSS, `RLIMIT_AS=10 GiB`, threads 2 — [`health/157_verify_all.log`](health/157_verify_all.log) |

**Scope distinction, stated explicitly.** `verify_all.sh` hashes the immutable
external item manifests (≈52 GB of `pert_data/`, `pert_model/`, `disease_data/`,
`drug_emb_model/`, `pheno_data/`). It does **not** hash any Stage 2 expression
payload, does not open the bases `.h5ad`/`.npy` files and is not a bases-product
integrity check. Stage 2 product integrity is covered by the metadata audit (A8,
manifest digests with expression payloads deliberately excluded) and, for the
accepted payload bytes, by the #156 E1 re-hash record — not by C3.

## D. Tests

13/13 individual files PASS from the branch worktree (46 tests, 1 documented
opt-in skip), threads pinned to 2, peak 233,960 KiB. Counts match the #156
record exactly. [`health/157_tests.log`](health/157_tests.log).

## E. Findings

| id | finding | disposition |
|---|---|---|
| F1 | `execution/build_identity.py` still defaulted `--native output/native-v1`, `--output output/identity-v1` and `--shared output/shared-v1` — all three retired paths, so the README's bare `python …/build_identity.py` could not run and would have written to a deleted path | **FIXED.** Defaults rebound to `output/gene-identity` (axis/crosswalk + product) and `log/historical/shared-v1-mappings` (the surviving historical v1 mappings the seal binds, per #156 D2 — *not* the current shared-expression mappings, so no re-seal). Proved with `--verify` against the accepted product: manifest `5f7517c3…` and genes `5877a3cd…` reproduce byte-identically, `FAM153B` index 5091, full-axis HGNC round trip PASS |
| F2 | the #156 audit's rejected-Tahoe check asserted the payload was *present*, which the removal necessarily breaks | **REWRITTEN** as a removal-record check: PASS only while the path is absent, the historical manifest still declares `6e3ac2ec…`, the `log/157-removals.json` entry agrees, and both rejection records survive. An unrecorded disappearance still fails. This is the recorded D6 exception, not a relaxed check |
| F3 | no empty old-path directory survived under `output/` | nothing to remove; recorded in `157-removals.json → empty_dir_scan` |
| F4 | `log/diagnostic/cohorts.csv` duplicates `output/tenv-hvg/cohorts.csv` | kept and reported; #157 is not authorized to drop accepted audit evidence |
| F5 | #156's G9 claim of "no output" from its old-path grep is not what the tree shows: the grep also matched `build_identity.py` (F1) and the retained `log/historical/shared-v1-mappings` name | recorded; G9 is re-evaluated here as A4 and now holds for real. #156's record keeps its bytes as history |

## F. Deviations

1. `check_paths.sh` was captured through `tail -40`; the per-literal
   total/resolved line was emitted above that window and is not reproduced. The
   ticket allowed exactly one run, so no second run was made to recover it. The
   PASS result, the 0-unresolved count and the whitelist counts are verbatim.
2. The documented `systemd-run` cap is unavailable (no user bus); the equivalent
   10 GiB `RLIMIT_AS` was applied by `execution/capped.sh`, as in #156.
3. The metadata audit resolves `log/157-removals.json` from the step root, falling
   back to the directory that holds the audit script. That lets a split checkout
   whose payload lives in the main checkout verify the removal record that
   travels with the code, while the record it writes stays free of worktree paths.

## G. Reused vs new evidence

Machine-readable ledger: [`157-evidence-ledger.json`](157-evidence-ledger.json).
