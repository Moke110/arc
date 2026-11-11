# #156 gates and acceptance evidence

Ticket #156, parent #154. Branch `codex/stage2-bases-156` based on `dev`
`22048777`. Worktree `/home/chang/arc/.scratch/worktrees/stage2-bases-156`; the
34 GB payload tree stayed in the main checkout `/home/chang/arc` and every tool
below took explicit `--step-root`/`--candidate`/`--figure`/`--records-dir`
paths. Nothing was copied into a worktree during the migration and nothing was
pushed or merged.

Follow-up acceptance fix: branch `codex/stage2-bases-156-fix` based on the merged
`dev` `9c3beba`. It removes the two superseded step-root placements, corrects the
audit's stable-set and duplicate-check semantics, drops the quoted fixture
literal from this record and fills the missing superseded-seal digest.

Environment: Python 3.13.13, matplotlib 3.10.9, pandas 2.3.3, NumPy 2.2.6;
R 4.5.3 / DESeq2 1.50.2 at `/home/chang/miniconda3/envs/biomni_e1/bin/Rscript`
(`4b08b95a…`). Heavy steps ran one at a time through `execution/capped.sh` with
`/usr/bin/time -v` and `OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS=2`.

**Cap substitution.** The documented `systemd-run --user --scope` form is not
reachable in this session (`Failed to connect to bus: No data available`).
`execution/capped.sh` detects that and applies the equivalent 10 GiB ceiling as
`RLIMIT_AS` locally; measured peaks were 22–853 MiB, so the ceiling was never
approached. See "Deviations".

## Bounded expression gates (E1–E5)

| id | check | outcome | measured | evidence |
|---|---|---|---|---|
| E1 | re-hash the accepted payloads once after the move | **PASS** | 1425 files, 20,880,197,658 B, 12.2 s wall, 22.4 MiB peak RSS | [`156-payload-digests.json`](156-payload-digests.json) |
| E2 | re-issue + verify the truth successor envelope, serial | **PASS** | create 243.4 s + verify 234.6 s = 7:59 wall, 878,732 KiB peak RSS | [`156-successor-reissue.json`](156-successor-reissue.json), [`records/156_truth_successor_build.json`](records/156_truth_successor_build.json), [`records/156_truth_successor_verify.json`](records/156_truth_successor_verify.json) |
| E3 | `verify_hvg_metrics.py` against the stable paths | **PASS** | 1:25 wall, 405,908 KiB peak RSS | [`records/156_e3_diagnostic_replay.json`](records/156_e3_diagnostic_replay.json) |
| E4 | render the panel from the CSV alone and compare PNG bytes | **PASS** | seconds, byte-identical `35d14266…` | [`156-figure-reproduction.json`](156-figure-reproduction.json) |
| E5 | `evaluate_predictions.py --successor … --diagnostic` end to end | **PASS** | 2.8 s, 365,728 KiB peak RSS | [`records/156_e5_evaluator_report.json`](records/156_e5_evaluator_report.json), [`records/156_e5_predictions_export.json`](records/156_e5_predictions_export.json) |
| E6 | refit only if E1–E5 contradict an accepted record | **not triggered** | — | no contradiction found |

E1 command:

```bash
execution/capped.sh /usr/bin/time -v python3 execution/verify_payload_digests.py \
    --step-root /home/chang/arc/stage2/01_bases --record …/log/156-payload-digests.json
```

It checks every accepted manifest's own seal, the two accepted diagnostic
evidence files, and every file those manifests declare (including both
`.h5ad` pairs, both shared `.npy` arrays and the truth statistic matrices). The
rejected historical Tahoe payload is explicitly skipped as a non-input.

E2: the accepted historical envelope's bytes (manifest `24b65aee…`) are
preserved at `log/historical/seals/truth-successor-envelope-v1/`. The published
envelope is the newly created and independently re-verified candidate, manifest
sha256 **`2798e3b8eb32027b221c18634bbb9e1dd4aeed46624786bd8099fc40fc790b6d`**
(record: [`156-successor-reissue.json`](156-successor-reissue.json); create and
verify outputs: `log/records/156_truth_successor_{build,verify}.json`). A first
candidate (`0e55fd98…`) was withdrawn before acceptance because it bound the
temporary ticket-worktree copy of the verification record; its bytes and the
defect are recorded in
[`records/156_rejected_successor_candidate.json`](records/156_rejected_successor_candidate.json).
The final envelope binds the canonical paths only:
`output/tenv-truth/manifest.json`, `output/historical-training-release/{manifest.json,tenv.h5ad}`,
`output/training-corpora/{manifest.json,tenv.h5ad}`, `output/tenv-raw-counts/tenv.h5ad`
(+ its producer manifest), `log/records/truth_complete_verify.json` and the
current `Rscript`. Create reported 1377 contrasts / 4206 samples; verify reported
`groups_compared: 1377`, `gene_count: 18164`. The re-issue record now records the
superseded envelope digest explicitly (`superseded_manifest_sha256`
`24b65aee…`), recovered from the preserved seal.

E3 (re-run against the final envelope, successor seal `2798e3b8…`): the replay recomputed the HVG source replay and the held-out isolation
intervention for all three folds, then exported the diagnostic panel to a fresh
candidate directory. Every plotted artifact was **byte-identical** to the
accepted evidence: `diagnostic_pr.png` `35d14266…`, `diagnostic_pr.csv`
`eaad870b…`, `diagnostic_pr.json` `fb364f05…`, `fold_{0,1,2}_groups.csv`
`7a543149…/dae8b0b7…/e635c406…`, and the per-fold macro AP values and eligible
group counts equal the accepted `acceptance.json`. The accepted evidence in
`log/diagnostic/` was not overwritten.

E4: the published CSV carries 7222 curve rows and 9 summary rows (three folds ×
three panels). The renderer was executed with `HOME`/`MPLCONFIGDIR` pointed at a
scratch directory and only `figure/data/diagnostic_tenv_pr.csv` as input; the
PNG digest equals the accepted `35d14266…`.

E5 (re-run against the final envelope): fold 0 measured-delta rows (1401 samples)
were exported in bounded blocks and scored through `SuccessorTruth` against the
re-issued envelope. Macro AP is
exactly the accepted fold 0 value (`0.07981072489453635`, positive
`0.3349457015461418`, negative `0.13737174931309615`) and the per-group CSV
digest is the accepted `7a543149…`.

## #156 publication gates

| id | outcome | evidence |
|---|---|---|
| G1 stable paths only | **PASS** | the migrated payload tree holds exactly the nine stable product directories plus the `output/api` code publication; the branch tracks `output/api/bases/**` and no file at an old path, and the two superseded step-root placements are gone. `find output -maxdepth 1 -type d` + `git ls-files`. |
| G2 no payload bytes changed | **PASS** | E1: 1425 files / 20.88 GB re-hashed to the declared digests. |
| G3 non-expression audit | **PASS** | [`156-metadata-binding-audit.json`](156-metadata-binding-audit.json), regenerated on the fix branch (split checkout, code root = branch tree): 44 binding + 21 structural checks, `non_pass: []`, one documented `EXPECTED_GAP` (the deleted shared-v1 seal that the historical release references). |
| G4 successor re-issued, not edited | **PASS** | E2 + preserved superseded bytes. |
| G5 evaluator through new paths only | **PASS** | E5; no code path falls back to an old path (see G9). |
| G6 one figure triple, CSV-only render | **PASS** | E4 + audit figure checks. |
| G7 step-root duplicates gone | **PASS** | the two superseded step-root placements are removed from the tracked tree (follow-up commit); the PNG lives at `figure/diagnostic_tenv_pr.png` (`35d14266…`, byte-identical to the preserved evidence copy `log/diagnostic/diagnostic_pr.png`), the raw CSV bytes at `log/diagnostic/diagnostic_pr.csv` (`eaad870b…`) and the extended plotted-values CSV at `figure/data/diagnostic_tenv_pr.csv`. |
| G8 evidence out of `output/` | **PASS** | `log/diagnostic/`, `log/records/`, `log/health/`, `log/historical/`; no `*.h5ad`/`*.npy` under `log/` (audit check). |
| G9 no runtime old-path reference | **PASS** | see the exact command and residual list below. |
| G10 README/AGENTS in the new vocabulary | **PASS** | `stage2/01_bases/README.md`, `AGENTS.md`, `execution/hvg_metrics.md`, `execution/truth_successor.md` reviewed; D1 placement documented. |
| G11 data inventory | **PASS** | `docs/data-inventory.md` lists every stable product with refreshed `du -sb` sizes, the api code and the figure/records roles. |
| G12 repo health | **PASS** | `check_tracked` PASS; `check_paths` PASS on the final branch tree (858 text files, 6254 literal occurrences, 0 unresolved); the migration-time run's two pre-existing hits are fixed and recorded (details below). `verify_all.sh` not run. |
| G13 existing science still binds | **PASS** | E1 digests + E3 byte-identical panel/AP + E5 identical macro AP and groups digest; no refit, no retrain. |
| G14 rejected Tahoe still present/rejected | **PASS** | audit check on `historical-training-release/tahoe_bulk.h5ad` digest `6e3ac2ec…`; no code path opens it. |

G9 command (run in the worktree), with the one-shot migration script and the
records that describe the move as the only expected hits:

```bash
git grep -nE 'identity-v1|native-v1|shared-v[12]|training-v[12]|truth-tenv-v[12]|tahoe-counts-v1|hvg-tenv-v1|hvg-metrics-diagnostic-v1|build_record' \
  -- stage2/01_bases ':!stage2/01_bases/log/records' ':!stage2/01_bases/log/health' \
     ':!stage2/01_bases/log/historical' ':!stage2/01_bases/log/155-*' ':!stage2/01_bases/log/156-*' \
     ':!stage2/01_bases/execution/migrate_products.py'
# no output
```

## Health checks

Re-verified from the follow-up branch `codex/stage2-bases-156-fix` (based on the
merged `dev`): the tracked tree now carries neither the step-root duplicates nor
a live-looking fixture literal in this record.

`script/check_tracked.sh` (fix worktree): **PASS** — 242 code files on
disk, 284 tracked, 0 untracked, 0 holes, 0.01 s, 4.4 MiB.
Record: [`health/156_check_tracked.log`](health/156_check_tracked.log).

`script/check_paths.sh` (fix branch tree; the main checkout owns the external
payloads, so the run used the documented temporary `GIT_INDEX_FILE` from
`git read-tree codex/stage2-bases-156-fix`, inside a private mount namespace
where the branch `stage2/` tree was bind-mounted read-only over the main
checkout's `stage2/`): **PASS** — 8–9 s wall, 6.3 MiB peak RSS,
858 text files scanned, 6254 literal occurrences / 248 distinct,
6224 resolved, 30 whitelisted (20 distinct), 0 unresolved.
Record: [`health/156_check_paths.log`](health/156_check_paths.log).

The migration run was red for two pre-existing hits, both reproduced on
untouched `dev`:

1. a fixture literal captured in the frozen
   `final_health/check_paths_corrected.stdout.log` — it never existed in the
   repository; that frozen file is exactly the evidence the whitelist exists
   for, so it stays whitelisted by path (`script/path_whitelist.txt`) and the
   match predicate was re-checked directly. This record previously quoted the
   literal itself; it now names the frozen file instead, so the repository path
   check no longer sees a live-looking path in a document.
2. `check_stage2_paths.sh` reported `unexpected Stage 2 entry: stage2/AGENTS.md`
   on pristine `dev` as well: its entry loop compared every `stage2/*` entry to
   `stage2/01_bases` without requiring a directory. The loop now skips
   non-directories (the step-contract document is not a step directory);
   `script/tests/stage2_paths_audit.test` and
   `script/tests/bases_scope_health.test` both PASS
   ([`health/156_health_script_tests.log`](health/156_health_script_tests.log)),
   and the standalone stage2 audit now PASSes
   ([`health/156_stage2_path_audit.log`](health/156_stage2_path_audit.log)).

## Tests

Each file was inspected and run individually (`python3 -m unittest discover -s
tests -p '<file>'`, threads pinned to 2); all PASS: `test_identity` (5),
`test_metrics` (10), `test_plots` (1), `test_raw_provenance` (1), `test_native`
(1), `test_hvg` (2), `test_corpora` (8), `test_significance` (3),
`test_evaluate_predictions` (1), `test_tahoe` (3), `test_truth_successor` (8),
`test_significance_worker` (1, skipped: opt-in `ARC_DESEQ2_TEST`),
`test_truth_successor_runtime` (2). The full suite was not run as one command.

## Reused vs new evidence

| claim | reused record | #156 re-ran it? |
|---|---|---|
| tenv raw counts 4510 wells, source replay | `native_tenv_verify.json` | no; digest re-hashed (E1) |
| Tahoe raw counts, 3388 shards / 95624334 cells | `tahoe_counts_source_acceptance.json` | no; payload digest re-hashed (E1) |
| shared expression zero-error replay | `shared_v2_verify.json` | no; payload digests re-hashed (E1) |
| training corpora build/verify and both source replays | `training_v2_{build,verify,tenv_replay,tahoe_replay}.json` | no; digests + E3 source replay |
| truth fit, 1377 contrasts, three R refits | `truth_complete_verify.json` | no refit; E2 replays it, E3 replays statistics |
| successor equivalence proof | old envelope `24b65aee…` | **yes** — re-issued and fully re-verified (E2) |
| HVG rankings and held-out isolation | `hvg_actual_{build,verify}.json` | **yes** — replayed (E3) |
| diagnostic panel values and visual QA | `diagnostic/{acceptance,visual_qa}.json` | **yes** — re-rendered (E3/E4); bytes identical |
| evaluator end-to-end | not previously recorded for the stable paths | **yes** (E5) |
| identity axis and mappings | `155-metadata-binding-audit.json` | new digest audit (G3) |
| layout, records, figure schema | #155 map and gates | **yes** — new audit (G3/G8) |

Machine-readable per-claim ledger: [`156-evidence-ledger.json`](156-evidence-ledger.json).

## Deviations and limitations

1. **Cap substitution**: `systemd-run --user --scope` was unreachable; the same
   10 GiB ceiling was applied as `RLIMIT_AS` by `execution/capped.sh`, and the
   measured peak RSS is reported for every heavy step. No job approached the cap.
2. **Split-checkout verification**: the tracked code/figure live in the ticket
   worktree, the payload in the main checkout; tools take explicit roots. G1 is
   therefore stated as two facts (payload directories in the data root, tracked
   code paths in the branch) rather than one `find` over a single checkout.
3. **`check_paths` re-run on the follow-up branch.** The migration run was not
   repeated before the merge; the follow-up branch drops the record's quoted
   fixture literal and runs the check against the final tracked tree (private
   bind mount of the branch tree; see Health checks). The plain pre-merge
   invocation against the main checkout still reads the merged `dev`
   working-tree copy of this record, which cannot be updated before the branch
   lands, so that invocation is expected to turn green at merge time.
4. The published `figure/data` CSV necessarily differs from the accepted
   evidence CSV (`eaad870b…`) because it now carries the displayed macro-AP
   values and the fixed-group definition; the original bytes are preserved and
   digest-checked at `log/diagnostic/diagnostic_pr.csv`.
5. **No stochastic tolerance is claimed**: all replays here are deterministic and
   reproduced bytes exactly; file integrity is reported as sha256, separate from
   metric values.
6. **Main-checkout state after the migration.** The main checkout must stay on
   `dev`, so its dev-tracked compact files at the old paths (140 paths: the
   step-root diagnostic pair, `execution/build_record/**` and the force-tracked
   files inside the moved products) were restored from `dev` HEAD after the
   migration to leave `git status` with no tracked changes. Those restored
   copies are `dev`'s state, not the branch's layout. The shipped branch tracks
   the new stable paths only and **deletes the two superseded step-root
   placements** (`diagnostic_tenv_pr.{png,csv}`); their bytes are preserved at
   `figure/diagnostic_tenv_pr.png` (`35d14266…`, also `log/diagnostic/diagnostic_pr.png`)
   and the raw plotted values at `log/diagnostic/diagnostic_pr.csv`
   (`eaad870b…`). The audit record was regenerated on the follow-up branch in
   split-checkout mode (payload root = main checkout, code root = branch tree);
   its duplicate check judges publishing on the code root, so `non_pass` is
   empty there. After the branch lands on `dev` the default single-checkout run
   (all roots = the main checkout) is expected to be identical, because the
   tracked copies no longer exist.
7. `log/records/truth_successor_acceptance.json` and
   `log/records/hvg_metrics_actual_acceptance.json` keep their original bytes and
   therefore name the *superseded* envelope digest (`24b65aee…`); the audit
   binds them to the preserved envelope bytes at
   `log/historical/seals/truth-successor-envelope-v1/` and the new digest is
   recorded in the re-issue record.
