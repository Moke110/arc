# Review-fix note — `stage2/tsg_view`

Code/test/evidence fixes for the two-axis review of `main...dev` (69 commits). No sealed
product byte was changed: `output/views`, `output/truth`, `output/report`, `figure/`, the
freezes, `formal_plan.json`, `formal_plan_manifest.json` and `run_records.json` are
byte-identical to `dev`. Only the code, tests, `log/build-03.json`, the new
`log/ticket-03-discarded-members.csv` and the regenerated (gitignored)
`output/runs/formal_verify.json` changed.

| # | Finding | Status | Where |
|---|---------|--------|-------|
| 1 | Full-arm anchor conflates seeds (spec c1) | fixed | `execution/formal.py:607-660`, `execution/tests/test_formal.py:FullArmAnchorSeedTest` |
| 2 | Ticket 03 item 4 discarded-member evidence missing (spec a3) | fixed | `execution/evidence.py`, `log/ticket-03-discarded-members.csv`, `execution/tests/test_evidence.py`, `log/build-03.json:727` |
| 3 | False 1.06 GB accumulator claim | fixed | `execution/views/__init__.py:688-699` |
| 4 | "Hermetic, no real payloads" claim | fixed (tests) | `execution/tests/{test_formal,test_matrix,test_learner_smoke,test_truth}.py`, `execution/tests/tsg_view_fixtures.py` |
| 5 | Baseline duplication smells | fixed, except `control_rule` (deferred) | `execution/common.py` + call sites |
| 6 | Evidence hygiene (spec c2) | this note | — |

## 1. Full-arm anchor seed

`matrix.FULL_ANCHORS` documents the sealed Step 04 anchors as equal-fold means at **seed
42**, but `verify` compared every `full` origin, seed 43 included, with
`matrix.compare_full_arm` and printed `PASS`. `formal_verify.json` therefore presented
seed-43 fits as reproducing sealed scores. `full_arm_evidence()` now splits the arm: only
seed-42 origins and fold means are compared with the anchors (`PASS`/`DRIFT`) and decide
the report `status`; every seed-43 origin and fold mean is recorded under
`full_arm.sensitivity` with `status: "sensitivity"`, `anchor_seed: 42`, the measured
value, the drift and a note that it claims no reproduction. `sensitivity_record()` refuses
to label the anchor seed itself as a sensitivity, so the split cannot silently invert.

## 2. Ticket 03 discarded-member evidence

The per-condition discarded count was only derivable from the payload's gitignored
`members.csv`. `execution/evidence.py` regenerates one row per condition of each
member-choice product from the published `members.csv` companions (never a payload
matrix), and the result is tracked as `log/ticket-03-discarded-members.csv` (55,531 rows:
solo tenv 1,377, solo tahoe 52,777, solo-deepest tenv 1,377). `test_evidence.py`
reconciles the committed table with the released population constants (tenv 1,320x3 /
48x4 / 9x6, Tahoe 43,036 singletons / 7,941x2 / 1,400x3 / 350x4 / 50x16; discarded 2,829
for each tenv arm and 12,541 for Tahoe), recomputes every count from the independently
built `merge` member table and proves regeneration is byte-identical. No manifest, payload,
freeze, truth product or run was rebuilt or edited.

## 3. Accumulator comment

float64 and int64 are both 8 bytes, so the single-pass `52,777 x 5,000` accumulator is
about 2.11 GB either way; the old comment described a gene-block two-pass the code does
not implement. The comment now states the real footprint and why float64 summation is
nevertheless exact (at most 16 small integer addends).

## 4. Test hermeticity

The tests no longer hard-code `/home/chang/arc`: the native source is repo-relative or
`TSG_NATIVE_SOURCE`, the gitignored payloads are repo-relative or `TSG_PAYLOAD_ROOT`, and a
class that genuinely needs one is `skipUnless`-guarded with a message naming the path and
the override. `test_truth`'s silent off-host skip is gone. In a bare worktree the suite is
192 tests with explicit skips; pointed at the materialised dependencies it is 192/192.

## 5. De-duplication

One implementation now lives in `execution/common.py` and is imported everywhere: the
canonical JSON `dumps`/`write_json`/`read_json`, the repo-root finder, `relative`,
`sha256`, the source/view order and per-view seed list, and the released-population
constants plus member-count histogram used by `views/merge.py`,
`views/full_pooled_label.py` and `truth.py`. The `_condition_identities` middle man, the
one-item tuple loop in `render.render`, the duplicate `learning_curves(runs)` call in
`report.figure_data` and `learner.TruthBlock`'s undocumented alternate columns are gone.
The emitted truth column names are unchanged (`drug_parent`, `members`).

### Deferred: `control_rule` token vocabulary

`CONTROL_RULE` is prose in four shapes ("mean of the matched_control_group vehicle rows,
each normalized ...", the merge normalize-of-sum sentence, the full-view-control sentence
for both member-choice arms and the pooled-label sentence). The review suggested one token
vocabulary (e.g. `mean-normalized` / `normalize-of-sum` / `view-control`) plus one
explanation string. That would change the `control_rule` value stored in every published
view manifest, and `output/views/**` is a sealed product whose manifest digests are pinned
by the freeze, the run configs and the acceptance receipt; rewriting the strings would
therefore break the byte-stability gate this pass is required to preserve, and it is a
semantic-payload change dressed as a refactor. Deferred to a round that may rebuild and
re-seal the views: it needs the manifests, freezes, runs and receipt regenerated together,
not a code-only edit.

## Code-identity consequence (must be carried into the final receipt)

The sealed plan binds two code-digest blocks, `code` and `learner_runtime`. The
de-duplication changes `matrix.py`, `learner.py` and `views/__init__.py`, so
`code_identity` in `formal_verify.json` now reports `CODE_DRIFT` with a section-tagged
`code_drift` list. `verify` records that drift instead of refusing to run (otherwise no
refactor could ever be verified); every non-code input change (freezes, truth, native
source, runtime environment, schedule) still fails closed, and `run`/resume keeps the
driver-file-only allowance. A *replay* of the sealed runs still refuses a moved learner
runtime (`learner.load`), which is the intended guard: the published run artifacts are
unchanged, but regenerating predictions would need the pre-refactor worker. The final
receipt pass must state this and regenerate the receipt against the new HEAD.

## Docs the final pass must apply

* `AGENTS.md:97` — replace "The tests are hermetic (synthetic fixtures, no payloads, GPU or
  network)." with a two-tier statement: the core suite is hermetic; the payload-dependent
  tests find the TranSiGen source and the gitignored payloads repo-relatively or through
  `TSG_NATIVE_SOURCE` / `TSG_PAYLOAD_ROOT`, and skip with a message naming the missing path.
  Give the full-coverage command with both variables set.
* `AGENTS.md`/`README.md` — mention `log/ticket-03-discarded-members.csv` (and
  `execution/evidence.py`) among the ticket-03 evidence; cite `log/REVIEW_FIXES.md` from the
  final report.
* If the docs repeat the 1.06 GB accumulator figure, replace it with the corrected bound.
