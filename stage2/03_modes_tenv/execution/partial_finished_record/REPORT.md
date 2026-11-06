# Finished-model partial evaluation, 2026-11-06

Explicit user amendment: “Amend the contract; score finished models now.” ADR 0042
permits this dated snapshot while the unchanged serial formal queue continues.
This is not final 63-run evaluation and never satisfies its completion barrier.
Scores must never affect training, checkpoint selection, hyperparameters, retries,
scheduling, admission or entrant changes. All final seal/evaluation/replay/delivery
and tracker closure obligations remain outstanding.

Captured at 2026-11-06T08:31:44.370120+00:00: **47 completed**; both references,
PRnet, CPA and chemCPA each 9/9; CMonge 2/9 (fold 0, seeds 17 and 42), descriptive
only; PerturbDiff 0/9 (one running). The other 15 identities were not dispatched.
Pending state is a read-only snapshot, not a fabricated terminal outcome.

The unchanged frozen `benchmark.evaluate` function supplies all group metrics,
well-to-group averaging, run means/sample SD, PR bands and fixed-order best-run
selection. An isolated function namespace substitutes only a validated partial
manifest gate, the read-only source-run root and rendering callback. Every
completed run matches frozen expected identity plus live original environment,
source, input, split, config, prediction and checkpoint hashes; validation repeats
after scoring. No shared module globals, source predictions or training files
are modified. Noncompleted identities are never scored or assigned zeros.
`predictions.json` deliberately has `sealed:false` and `final_barrier_passed:false`.
The regression tests demonstrate that the canonical barrier rejects it.

## Actual results

Nine-run entries are the primary partial comparison. The two CMonge runs share
fold 0, so they cannot establish a whole-fold or complete-entry ranking. The
best-run display is post-hoc test-selected among available completed runs only,
can select different folds, and does not replace the primary nine-run summary.
Deterministic reference seed copies are not independent fits. Undefined no-change
correlations remain NaN with zero metric-defined counts.

| Entry | Completed | AP positive | AP negative | PCC | RMSE |
|---|---:|---:|---:|---:|---:|
| no_change | 9/9 | 0.0055054 | 0.0055054 | nan | 0.14445 |
| delta_mean | 9/9 | 0.189543 | 0.17367 | 0.229553 | 0.140614 |
| prnet | 9/9 | 0.182874 | 0.0970652 | 0.257316 | 0.183764 |
| cpa | 9/9 | 0.0711623 | 0.0722207 | 0.293605 | 0.142937 |
| chemcpa | 9/9 | 0.0429868 | 0.047719 | 0.183627 | 0.178152 |
| cmonge | 2/9 | 0.208161 | 0.192108 | 0.22563 | 0.141717 |
| perturbdiff | 0/9 | nan | nan | nan | nan |

## Artifacts and execution

Visible figures: `output/metrics.{png,svg}`, `output/pr_panels.{png,svg}` and
`output/best_run.{png,svg}`. **Final full-terminal delivery must overwrite these
same six visible paths.** The dated snapshot remains separately under
`output/partial_finished/20261106_snapshot/`; all five metric CSVs and the curve
NPZ are local there. This record retains the compact snapshot manifest, checked
evaluation/source hashes, publication hashes and summary/run/band/selection CSVs.
`output/partial_publication.json` records the visible figure provenance and
explicit overwrite policy. No canonical final seal/evaluation was written.

Reproduction: `command.sh` (fresh snapshot directory required; it intentionally
refuses overwriting a snapshot). Input verification and completed-run hashing
precede scoring. Workers=2; all BLAS/OpenMP threads=2; systemd MemoryMax=10G and
MemorySwapMax=8G. Preflight available RAM=38GiB. This was the second heavy CPU
job beside the unchanged GPU worker; it has ended and released its slot.
`evaluation.log`: exit 0, 2:49.58 wall time, **477680 KiB peak RSS**.
`tests.log`: two artifact-boundary tests pass, 4.87s, **207720 KiB peak RSS**.
`delivery.log` records lightweight final rendering/publication resource use.
No GPU work, inference, refitting, scheduler changes or score-driven decisions.

Visual QA: inspected all three real PNGs. Eight metric panels show finite-run
means, dots/sample SD, per-metric defined counts, undefined values, orange
incomplete region and deterministic-copy note. PR panels use red/blue directions
with min/max bands and explicitly tinted incomplete panels. Best-run legend
shows selected fold/seed, all completion counts, dashed incomplete curves and
available-run test-selection disclosure. All figures visibly show PARTIAL,
snapshot date/counts and primary/descriptive scope. SVG text remains editable
text; publication hashes match all six visible copies. Focused tests also reject
tampered scores, source predictions, frozen identities and invalid publish paths.

Independent narrow Standards/Spec review: **PASS**, no actionable blocker.
Reviewer independently checked exact frozen evaluator reuse, read-only source and
identity boundaries, all scored-artifact hashes, six visible publication hashes,
nonsealed 47/63 snapshot and final-barrier-false status, and visually inspected
all three real PNGs. Initial presentation finding (incomplete-entry separation
and available-run best-selection disclosure) was resolved before commit.
Generated CSV CRLF and Matplotlib SVG whitespace are preserved byte-for-byte
with their provenance hashes; source-only whitespace checks pass.
