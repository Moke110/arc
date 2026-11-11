# #156 decisions

Ticket #156, parent #154. Decisions below are the coordinator's and are recorded
here with their rationale and the mechanical evidence; #156 does not re-open them.
Authorities: ADR 0062, the `stage2/AGENTS.md` step contract and
[`155-migration-map.json`](155-migration-map.json).

## D1 — option A, contract-literal: one public implementation under `output/api/`

**Decision.** The single public implementation is published at
`stage2/01_bases/output/api/bases/` (the former `execution/bases/` package,
moved with `git mv`, bytes unchanged). Every build, audit, figure-provenance and
test command imports it through `execution/public_api.py`, which puts
`output/api` on `sys.path`. There is exactly one *live* copy of each module: no
module is duplicated as a second implementation and no second import root exists.
The accepted successor envelope carries a byte-identical, digest-sealed provenance
copy of `deseq2_contrasts.R` at `output/tenv-truth-successor/`, which
`truth_successor.py` re-hashes to prove the worker that produced the envelope;
that sealed copy is load-bearing evidence, not a duplicate, so keep it.

**Why A.** `stage2/AGENTS.md` states that `output/` holds public metric/reader
code and that *public implementations live here and execution imports them*;
ADR 0062 repeats "public readers and metric implementations belong to the output
contract". Option B (descriptor under `output/` pointing at `execution/`) only
partly satisfies that sentence and still forces readers to add `execution/` to
`sys.path`. Option C needed an ADR edit and is not #156's decision to take.

**Mechanics.** `.gitignore` ignores the *contents* of `output/` rather than the
directory itself and re-includes `output/api/**`, so the code is tracked while
34 GB of payload stays local. `script/check_tracked.sh` therefore still has zero
holes, and a future module under `output/api` cannot be silently swallowed. The
R DESeq2 worker moved with the package; `truth_build.py` seals
`sha256(Path(__file__))` and `sha256(significance.py)`, and both digests are
unchanged because this was a rename.

> **Correction (2026-11-12, #157 code review S1).** The Mechanics paragraph above
> is wrong about the mechanism: it says `.gitignore` "ignores the *contents* of
> `output/` rather than the directory itself" and that a future module under
> `output/api` "cannot be silently swallowed". The pre-existing `output/` rule
> excludes the `stage2/01_bases/output` directory itself, and git cannot re-include
> a path inside an excluded directory, so `git add` rejected **every** new file
> under `output/api` (reproduced: `git check-ignore -v` named `.gitignore:62:
> output/` and `git add -n` refused it). The decision stands; the mechanics are
> fixed in `.gitignore` by re-including the directory (`!…/output/`) *before*
> excluding its contents and re-including the api subtree. Verified with a real
> probe: `git add -n` accepts a new `output/api/bases/newmod.py`, and
> `git check-ignore -q` still reports `output/tenv-truth/log2_fold_change.npy` as
> ignored. `script/check_tracked.sh` proves the current tree only, not future
> adds, so this probe — not that script — is the regression check for the
> mechanism.

**Consequence for split checkouts.** During #156 the tracked code/figure live in
the ticket worktree while the 34 GB payload lives only in the main checkout.
Tools take explicit `--step-root`/`--figure`/`--code-root`/`--output` paths
instead of copying data; no payload was copied into a worktree.

## D2 — do not re-seal the identity product

**Decision.** `identity-v1` keeps its sealed bytes (now at
`output/gene-identity/`); the historical shared-v1 mapping CSVs are preserved at
`log/historical/shared-v1-mappings/`; the old-to-new binding is recorded in the
migration record and in `output/api/bases/identity.py`'s callers. The deferred
option of refreshing the 229 Tahoe `native_symbol` strings against the current
shared-v2 mappings is **not** taken.

**Why.** Re-sealing would change a public evidence interface (229 evidence
strings and the manifest digest) with no axis change, and #155 recorded it as a
user decision, not an implementation tidy-up. The identity manifest binds input
digests only, so moving the files leaves every declared digest valid (map
principle P2).

## D3 — no deletion of accepted, rejected or superseded payloads

**Decision.** #156 moves; #157 deletes. Every payload that the map calls
disposable was **moved** out of `output/`, not removed:
`output/truth-tenv-v1.building/` → `execution/checkpoints/truth-tenv-building/`;
`tahoe-raw-counts/.checkpoints/` → `execution/checkpoints/tahoe-raw-counts/`.
The rejected Tahoe training payload stays in place inside
`historical-training-release/` and stays marked rejected. The superseded truth
successor envelope's bytes are preserved at
`log/historical/seals/truth-successor-envelope-v1/`.

## Interpretations recorded for review (not re-opened coordinator decisions)

1. **Step-root duplicates (gate G7).** `diagnostic_tenv_pr.png` was moved to
   `figure/diagnostic_tenv_pr.png` and `diagnostic_tenv_pr.csv` to
   `figure/data/diagnostic_tenv_pr.csv`. The original bytes are preserved as
   audit evidence: the PNG is byte-identical (`35d14266…`) at
   `log/diagnostic/diagnostic_pr.png`, and the *original* plotted CSV bytes
   (`eaad870b…`) at `log/diagnostic/diagnostic_pr.csv`. The published
   `figure/data` CSV extends those bytes with the displayed macro-AP values
   (finding F6), so its digest necessarily differs from the accepted evidence
   digest; the evidence copy keeps that digest checkable. The two superseded
   step-root placements were still tracked by the merge commit; the follow-up
   fix (`codex/stage2-bases-156-fix`) removes them from the tracked tree, which
   is the migration completing, not a payload deletion.
2. **`output/api` joins the stable set (D1).** Gate G1 was written before D1 was
   decided; the stable set is the map's products **plus** the code directory
   `output/api/`. G1 is evaluated as: the migrated payload tree contains exactly
   the map's stable product directories, and the branch tree tracks
   `output/api/bases/**` with no tracked file at an old path.
3. **Old-path literals (gate G9).** Two classes of file necessarily keep old path
   strings: the one-shot migration table in `execution/migrate_products.py`
   (its *inputs*, not a runtime default) and the #156 evidence records that
   describe the move. No runtime default, import or CLI in `execution/`,
   `output/api/` or `figure/` names an old product path; the gate command in
   [`156-gates.md`](156-gates.md) prints the residual list and its size (exactly
   those files).
4. **Tahoe mapping duplication removed without deletion.** The diagnostic
   directory's `cohorts.csv` was a byte-identical duplicate of the HVG product
   copy. It is retained as audit evidence at `log/diagnostic/cohorts.csv`; the
   single *product* copy remains `output/tenv-hvg/cohorts.csv` (#157 may drop the
   evidence duplicate).
