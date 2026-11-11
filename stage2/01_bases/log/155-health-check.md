# #155 health-check record

Exact scope of every health check run for this ticket, with the commands, the
working tree, and the honest outcome. Nothing here claims a check that was not
run or a verdict that was not produced.

Worktree: `/home/chang/arc/.scratch/worktrees/stage2-bases-155` (branch
`codex/stage2-bases-155`, based on `dev` `a53c9ff`). All #155 files are additions
under `stage2/01_bases/log/`; no product, script or tracked evidence file was
modified.

## check_tracked.sh — PASS

```bash
cd /home/chang/arc/.scratch/worktrees/stage2-bases-155
/usr/bin/time -f "wall=%es rss=%MKiB" script/check_tracked.sh
```

```
scan: 232 code files on disk (excluding virtualenvs, source/, __pycache__/, .git/, workspace archives); 274 tracked by git
untracked: 0  (documented expected-ignored: 0, holes: 0)
RESULT: PASS — all scanned executable files are tracked (or documented expected-ignored)
exit: 0        wall=0.01s   peak RSS=4,420 KiB
```

Scope: every `*.py`, `*.sh` and `*.R` in this worktree except virtualenvs,
`source/`, `local/`, `__pycache__/`, `.git/` and workspace archives. It reads
file names and `git ls-files`, never file contents, so it cannot hash an
expression payload. The one new executable file, `log/155-reproduce.py`, is
tracked, so the new `log/` directory introduces no hole and no ignore-rule
exception.

## check_paths.sh — NO VERDICT CLAIMED

Two bounded attempts in this worktree, both killed by the cap:

```bash
timeout -s KILL 90  bash -x script/check_paths.sh   # killed, 137
timeout            script/check_paths.sh            # background, killed at 240 s
```

Reasons, measured rather than guessed:

1. **Cost.** The script extracts 724 `pert_data/…`, `pert_model/…`,
   `disease_data/…`, `drug_emb_model/…` literals from tracked files and, for
   each occurrence and each of the 53 `script/path_whitelist.txt` lines, runs two
   subprocesses (`grep`, `sed`) inside `is_whitelisted()`. That is roughly 77 000
   process spawns; the run did not finish in 90 s or in a further 240 s. The
   reused health record's 9.99 s figure was measured at `b17049c`, before the
   fail-closed literal/whitelist revisions at `31b1560`/`a53c9ff`, so it is not a
   comparable cost.
2. **Environment.** An isolated worktree has no external `source/` payloads:
   `pert_data/Annotation/` and the other item directories are present only as
   tracked scaffolding. Any external literal therefore cannot resolve there even
   if the scan completes. This is by design, not a #155 regression.

Targeted predicate for the files #155 adds (same rule the script applies — literal
must resolve against the repository root), evaluated in the checkout that owns
the external data:

| literal in a #155 record | file | resolves at `/home/chang/arc` |
|---|---|---|
| `pert_data/Annotation/source/hgnc_complete_set.txt` | `155-inventory.md`, `155-reproduce.py` | EXISTS |
| `pert_data/GSE262419/source/` | `155-inventory.md` | EXISTS |
| `pert_data/Tahoe-100M/source/` | `155-inventory.md` | EXISTS |
| `pert_data/Tahoe-100M/source/metadata/gene_metadata.parquet` | `155-evidence-ledger.json` | EXISTS |

No other external-category literal appears in the new records.

## Reused health evidence (not re-run by #155)

`execution/build_record/final_health/summary.json` and
`final_health/evidence_manifest.json` record, at commits `b17049c`/`aaa48b5`:

- `check_paths`: PASS, 758 files scanned, 0 stderr bytes, 9.99 s wall.
- `check_tracked`: PASS, 244 files scanned, 0 stderr bytes.
- `verify_all`: exit 1 with four runtime item errors, retained exactly; the
  invalid-prior-verdict ledger keeps the earlier false passes visible.

#155 re-read those records and re-checked their internal digest bindings in
[`155-metadata-binding-audit.json`](155-metadata-binding-audit.json); it did not
re-run the checks and does not extend their verdict to the new layout.

## Guidance for #156

- Run `check_paths.sh` once, in the main checkout (where the external payloads
  exist), with a generous timeout, and record wall time and peak RSS. Expect the
  subprocess-bound cost above; do not run it twice in one session.
- Run `check_tracked.sh` in each worktree; it is cheap and layout-local.
- Do not run `script/verify_all.sh`. It hashes ~52 GB of external data, it was
  already recorded as failing for unrelated item-runtime reasons, and ADR 0062's
  audit does not need it.
