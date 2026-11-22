# tsg_view — does within-condition replicate structure explain the tenv–Tahoe gap?

A non-numbered **study step** for
[ADR 0082](../../docs/adr/0082-view-variants-tenv-tahoe-replicate-structure.md), routed from
the [Stage 2 step contract](../AGENTS.md#step-contract) and using the language of the
[glossary](../../GLOSSARY.md). It builds expression-view variants from the Bases integer
corpus under the ADR 0077 axis-local recipe, freezes each variant into the Step 04 input
format, fits native TranSiGen under the sealed `no-dose` cold-drug-dose folds, and scores
every variant on one frozen, view-independent truth. The question is how much of the
tenv–Tahoe Δ gap within-condition replicate structure in the labels or profiles explains.
The step publishes comparison evidence and **no model weights**, and it changes no product
of Steps 01–06.

## Pre-registration

Support for the hypothesis is a tenv Δ PCC gain from `merge` or `full-pooled-label` of at
most about 0.12 that stays below the single-sample ceiling near 0.38, with the tenv–Tahoe
gap still present. A larger gain, or one crossing the ceiling, revises the ceiling estimate.
Tahoe is the manipulation control because 81.5% of its conditions are unchanged by pooling
or drawing.

This block is the pre-registration recorded before the matrix was scored; it is pinned
verbatim as [log/acceptance/pre_registration.md](log/acceptance/pre_registration.md), and
its sha256 is `pre_registration_sha256` in the sealed
[acceptance receipt](log/acceptance/receipt.json). The verdict below is written against it.

## Results

The matrix is 45/45 origins complete: seven main views × three cold-drug-dose folds × seeds
42/43 plus three tenv `solo-deepest` fits at seed 42. `formal.py verify` reported
`status: PASS` with plan sha256 `2f59a870…`; the sealed receipt confirms 45 complete / 0
pending origins, exact replays and a CSV-only figure re-render. The headline estimand is the
**equal-fold mean** (mean of the three folds, then across seeds where reported); the primary
truth is per-profile Δ. Fold spread and seed spread are separate quantities and are never
pooled into one error bar. The tenv suite was published first; the Tahoe panels and the
cross-source `views_effect` were added afterwards from the same figure CSVs by a
deterministic renderer, and both passes are recorded in
[log/build-07.json](log/build-07.json).

### Primary truth — tenv suite (headline table, equal-fold mean)

| view | seed | ΔPCC | ΔSpearman | ΔRMSE | ΔR² |
|---|---|---|---|---|---|
| `full` | 42 | 0.2655 | 0.2861 | 0.4213 | -0.0824 |
| `full` | 43 | 0.2645 | 0.2899 | 0.4226 | -0.1021 |
| `merge` | 42 | 0.2324 | 0.2617 | 0.4367 | -0.1953 |
| `merge` | 43 | 0.2349 | 0.2649 | 0.4373 | -0.1902 |
| `solo` | 42 | 0.2074 | 0.2537 | 0.4543 | -0.3711 |
| `solo` | 43 | 0.2060 | 0.2573 | 0.4559 | -0.3942 |
| `full-pooled-label` | 42 | 0.2104 | 0.2596 | 0.4474 | -0.3017 |
| `full-pooled-label` | 43 | 0.2038 | 0.2556 | 0.4493 | -0.3130 |
| `solo-deepest` (sensitivity) | 42 | 0.2556 | 0.2944 | 0.4578 | -0.3343 |

Equal-fold mean over three folds per seed; `solo-deepest` is a single-seed sensitivity arm.
Source: [metrics_primary.csv](output/report/metrics_primary.csv), with fold and seed spread
in [fold_spread.csv](output/report/fold_spread.csv) and
[seed_spread.csv](output/report/seed_spread.csv).

### Primary truth — Tahoe suite (manipulation control, equal-fold mean)

| view | seed | ΔPCC | ΔSpearman | ΔRMSE | ΔR² |
|---|---|---|---|---|---|
| `full` | 42 | 0.4226 | 0.3584 | 0.1548 | -0.2811 |
| `full` | 43 | 0.4217 | 0.3569 | 0.1545 | -0.2789 |
| `merge` | 42 | 0.4093 | 0.3484 | 0.1571 | -0.3280 |
| `merge` | 43 | 0.4096 | 0.3471 | 0.1568 | -0.3243 |
| `solo` | 42 | 0.4142 | 0.3521 | 0.1569 | -0.3259 |
| `solo` | 43 | 0.4154 | 0.3509 | 0.1560 | -0.3061 |

Tahoe has no `full-pooled-label` arm: it is the manipulation control because 81.5% of its
conditions are single-profile, so pooling or drawing changes almost nothing by construction.

### Cross-source effect (`views_effect`)

Paired against `full` at seed 42, the tenv–Tahoe Δ PCC gap is **0.1572** for `full`,
**0.1769** for `merge` and **0.2068** for `solo`: every manipulation *widens* the gap rather
than closing it. [views_effect.png](figure/views_effect.png) shows the paired `view − full`
differences with a fold+seed interval, with `solo-deepest` in its own panel. Tables:
[paired_vs_full.csv](output/report/paired_vs_full.csv),
[paired_vs_full_solo_deepest.csv](output/report/paired_vs_full_solo_deepest.csv).

### Condition-mean truth

The same runs scored against the replicated-condition mean truth, ΔPCC per view × seed:

| view | tenv s42 | tenv s43 | Tahoe s42 | Tahoe s43 |
|---|---|---|---|---|
| `full` | 0.3599 | 0.3584 | 0.4179 | 0.4166 |
| `merge` | 0.3126 | 0.3109 | 0.3983 | 0.3972 |
| `solo` | 0.2995 | 0.3004 | 0.3992 | 0.4059 |
| `full-pooled-label` | 0.3505 | 0.3456 | — | — |
| `solo-deepest` | 0.3086 | — | — | — |

Test condition counts are tenv **153 / 153 / 153** and Tahoe **595 / 700 / 1,800** for
folds 0 / 1 / 2. The replicated-condition population is tenv 1,377 conditions / 4,206
profiles (3: 1,320, 4: 48, 6: 9) and Tahoe 9,741 / 22,282 (2: 7,941, 3: 1,400, 4: 350,
16: 50). **Caveat:** that population is *not* a random subset of its source — repeats
concentrate on multi-plate, multi-article conditions — so the condition-mean rows answer a
different question from the primary truth and are not a substitute for it. Source:
[metrics_condition_mean.csv](output/report/metrics_condition_mean.csv).

### Per-dose breakdown

Each fold holds exactly one dose, so the per-dose panel is that fold's sample-equal score,
not a dose average. **Dose and fold are confounded by design** and the two effects cannot be
disentangled here. ΔPCC per held dose (seed 42 / seed 43):

| view | tenv 0.2 µM (fold 0) | tenv 1 µM (fold 1) | tenv 10 µM (fold 2) |
|---|---|---|---|
| `full` | 0.2731 / 0.2787 | 0.2778 / 0.2754 | 0.2454 / 0.2394 |
| `merge` | 0.2529 / 0.2501 | 0.2278 / 0.2355 | 0.2166 / 0.2191 |
| `solo` | 0.2007 / 0.2137 | 0.2387 / 0.2270 | 0.1828 / 0.1772 |
| `full-pooled-label` | 0.2180 / 0.2060 | 0.2240 / 0.2210 | 0.1892 / 0.1845 |
| `solo-deepest` | 0.2855 | 0.2424 | 0.2389 |

| view | Tahoe 0.05 µM (fold 0) | Tahoe 0.5 µM (fold 1) | Tahoe 5 µM (fold 2) |
|---|---|---|---|
| `full` | 0.4487 / 0.4468 | 0.3981 / 0.3916 | 0.4210 / 0.4267 |
| `merge` | 0.4201 / 0.4206 | 0.3816 / 0.3831 | 0.4262 / 0.4251 |
| `solo` | 0.4242 / 0.4214 | 0.3881 / 0.3921 | 0.4303 / 0.4327 |

Source: [metrics_per_dose.csv](output/report/metrics_per_dose.csv) (all four metrics);
panels: [metrics_per_dose_tenv.png](figure/metrics_per_dose_tenv.png),
[metrics_per_dose_tahoe_bulk.png](figure/metrics_per_dose_tahoe_bulk.png).

### Each score beside its Δ noise floor and single-sample ceiling

The noise floor is the population-SD median of Δ over replicated conditions (headline); the
single-sample Δ PCC ceiling is the profile-weighted mean of `sqrt(rho_c)`, explicitly a
**heuristic diagnostic, not a measurement-error model**, and the input conditions are not a
random subset of the source.

| view | ΔPCC | ceiling | ΔPCC − ceiling | ΔRMSE | noise floor | ΔRMSE − floor |
|---|---|---|---|---|---|---|
| tenv `full` | 0.2650 | 0.3803 | -0.1153 | 0.4220 | 0.3426 | +0.0794 |
| tenv `merge` | 0.2337 | 0.3803 | -0.1466 | 0.4370 | 0.3426 | +0.0945 |
| tenv `solo` | 0.2067 | 0.3803 | -0.1736 | 0.4551 | 0.3426 | +0.1125 |
| tenv `full-pooled-label` | 0.2071 | 0.3803 | -0.1731 | 0.4483 | 0.3426 | +0.1058 |
| tenv `solo-deepest` | 0.2556 | 0.3803 | -0.1247 | 0.4578 | 0.3426 | +0.1153 |
| Tahoe `full` | 0.4222 | 0.4283 | -0.0061 | 0.1547 | 0.0703 | +0.0844 |
| Tahoe `merge` | 0.4095 | 0.4283 | -0.0188 | 0.1570 | 0.0703 | +0.0867 |
| Tahoe `solo` | 0.4148 | 0.4283 | -0.0135 | 0.1565 | 0.0703 | +0.0862 |

Tenv `full` sits **0.1153** below the ceiling. Tahoe `full` sits **0.0061** below its own
heuristic estimate — near that estimate, *not* near the ADR 0082 band: the estimate itself
(0.4283) is 0.0017 below the band floor of 0.43, so this is a moved heuristic estimate,
**not agreement with the 0.43–0.52 band**, and it was not tuned. The sample-SD sensitivity
variants are 0.4193 (tenv) / 0.0952 (Tahoe) for the noise floor. Sources:
[diagnostics_beside_scores.csv](output/report/diagnostics_beside_scores.csv),
`output/truth/<source>/diagnostics.json`. Panels:
[metrics_tenv.png](figure/metrics_tenv.png),
[metrics_tahoe_bulk.png](figure/metrics_tahoe_bulk.png).

### Paired view-vs-`full` differences

Paired differences over folds and seeds (mean and 95 % interval assembled from the fold
standard error and the seed standard error; fold and seed spread are never pooled into one
bar):

| view | ΔPCC | ΔSpearman | ΔRMSE | ΔR² |
|---|---|---|---|---|
| `merge` | -0.0313 [-0.0459, -0.0167] | -0.0247 [-0.0383, -0.0111] | +0.0151 [+0.0131, +0.0170] | -0.1005 [-0.1364, -0.0645] |
| `solo` | -0.0583 [-0.0731, -0.0436] | -0.0325 [-0.0392, -0.0258] | +0.0331 [+0.0235, +0.0427] | -0.2904 [-0.3789, -0.2018] |
| `full-pooled-label` | -0.0579 [-0.0665, -0.0492] | -0.0304 [-0.0437, -0.0170] | +0.0264 [+0.0216, +0.0312] | -0.2151 [-0.2752, -0.1549] |

`full-pooled-label` exists only for tenv. **`solo-deepest` is reported separately** as a
single-seed (42) three-fold sensitivity: ΔPCC −0.0099 [−0.0372, +0.0174], ΔSpearman +0.0083
[−0.0174, +0.0340], ΔRMSE +0.0365 [+0.0275, +0.0454], ΔR² −0.2519 [−0.3363, −0.1674]. Its
ΔPCC interval includes zero and it is not evidence of a gain. On Tahoe the same paired
differences are small: `merge` ΔPCC −0.0127 [−0.0296, +0.0042], `solo` −0.0074
[−0.0266, +0.0118].

### `Mean_cell` baseline

Step 04 axis-local `Mean_cell` scored on the same frozen cold-drug-dose rows (provenance and
receipt pinned in [report_manifest.json](output/report/report_manifest.json)):

| source | ΔPCC | ΔSpearman | ΔRMSE | ΔR² |
|---|---|---|---|---|
| tenv | 0.2270 | 0.2490 | 0.4418 | -0.7838 |
| tahoe_bulk | 0.4672 | 0.4170 | 0.1411 | -0.0045 |

Tenv `full` (0.2650 / 0.4220) beats `Mean_cell` on ΔPCC and ΔRMSE; Tahoe `full`
(0.4222 / 0.1547) is below `Mean_cell` on ΔPCC and close on ΔRMSE. Table:
[baseline_mean_cell.csv](output/report/baseline_mean_cell.csv).

### The `full` arm reproduces the sealed `no-dose` anchor

The `full` view is a digest-equal alias of the accepted 5k release, so its arm must
reproduce the sealed Step 04 `no-dose` cold-drug-dose scores. It does, to floating-point
noise, well inside the declared tolerance (`fold_mean` |Δ PCC| ≤ 0.01, |Δ RMSE| ≤ 0.005):

| source | sealed ΔPCC anchor | measured | drift | sealed ΔRMSE anchor | measured | drift |
|---|---|---|---|---|---|---|
| tenv | 0.2654588998492864 | 0.26545889984107024 | -8.22e-12 | 0.421337591468443 | 0.4213375914365656 | -3.19e-11 |
| tahoe_bulk | 0.4226328493964055 | 0.42263284939455525 | -1.85e-12 | 0.15479625072122077 | 0.1547962507163815 | -4.84e-12 |

The anchors come from `stage2/04_tsg_dose/execution/historical/report_axis_local/data/round1_summary.csv`;
the comparison is pinned in the receipt's `formal_verify` block.

The anchor comparison is **seed 42 only**. The seed-43 `full` rows are recorded in the same
`formal_verify.json` under `full_arm.sensitivity` with `status: "sensitivity"` and
`anchor_seed: 42`; they are reported as fold/seed spread, claim no anchor reproduction and
never decide the report `status`.

## Verdict against the pre-registration

**Not supported.** No manipulation produced a positive tenv Δ PCC gain: `merge` −0.0313,
`full-pooled-label` −0.0579 and `solo` −0.0583 (paired equal-seed means; every fold+seed
interval lies below zero), and `solo-deepest` −0.0099 at seed 42. The single-sample ceiling
was not crossed — every tenv arm remains below 0.3803 — and the tenv–Tahoe Δ PCC gap widened
from 0.1572 (`full`) to 0.1769 (`merge`) and 0.2068 (`solo`). Pooling labels and pooling
profiles therefore do not capture tenv's headroom under this protocol, and within-condition
replicate structure in the labels or profiles is not an explanation for the tenv deficit.
Tahoe moved only −0.0127 / −0.0074, consistent with its 81.5 % single-profile conditions.
Because the observed gain is not in the pre-registered direction, the ceiling estimate is
not revised by this outcome.

What the manipulation changed: tenv Δ PCC fell, Δ RMSE rose and Δ R² fell under every arm
(Δ PCC −0.031 to −0.058; Δ RMSE +0.015 to +0.033; Δ R² −0.10 to −0.29); Tahoe barely moved.
What it did not change: the view-independent truth, the Δ noise floor, tenv's distance below
the ceiling (0.1153 for `full`) and the fact that Tahoe sits near its own ceiling (0.0061).

## Post-review changes and code identity

The code, tests and ticket-03 evidence were revised after the two-axis review; the changes,
finding-by-finding status and the one deferral are recorded in
[log/REVIEW_FIXES.md](log/REVIEW_FIXES.md). No sealed product byte changed: `output/views`,
`output/truth`, `output/report`, `figure/`, the freezes, `formal_plan.json` and
`run_records.json` are byte-identical, and only the regenerated (gitignored)
`output/runs/formal_verify.json` moved. Because the refactor touches `formal.py`,
`learner.py`, `matrix.py` and `views/__init__.py`, that record reports
`code_identity.status = CODE_DRIFT` with six section-tagged moved digests: `verify` records a
code-only move instead of refusing it, while `run`/resume keeps the driver-file-only
allowance, and any non-code input change still fails closed. A fresh-process replay refuses a
moved learner runtime by design (regenerating predictions would need the pre-refactor
worker), so the sealed runs' replay evidence is the in-run, same-runtime evidence recorded in
`output/runs/formal_run_report.json` and the per-origin run records. A filtered `verify` that
selects no seed-42 `full` arm has no anchor comparison to make, so it now fails closed instead
of reporting a vacuous `PASS`; `full_arm.comparisons` records the count. The acceptance receipt
is regenerated after the review fixes; the re-pin is recorded in
[log/build-08.json](log/build-08.json) under `final_refresh`.

## Limitations and confounders not attributed

- **Control recipe.** `merge` pools the union of member control groups by summing retained
  integer counts, a declared generation that differs from ADR 0077 (a singleton condition
  keeps its X but its control deepens, so its Δ changes). It is never mixed into the `full`
  arm, but it is one candidate explanation for part of `merge`'s loss.
- **Dose and fold are confounded by design** (one held dose per fold); no dose effect is
  claimed.
- **Cell-line diversity, corpus size, assay and control structure differ together** between
  tenv and Tahoe and are not separated by this design.
- **The ceiling is a heuristic**, not a measurement-error model; a measurement-error model
  could replace it later and would revise the ceiling. Tahoe's heuristic estimate (0.4283)
  already sits 0.0017 below ADR 0082's 0.43–0.52 band; that is a reported deviation of the
  estimate, not agreement with the band.
- **The condition-mean truth's replicated-condition population is not a random subset** of
  its source; repeats concentrate on multi-plate, multi-article conditions.
- **Seed 43 is available for every main arm**, but `solo-deepest` is single-seed by design;
  seed spread is reported separately and never merged into fold error bars.
- **Tahoe is a manipulation control, not a replication**: 81.5 % of its conditions are
  single-profile, so its near-flat response confirms the manipulation reached the intended
  population but cannot test the hypothesis on its own.

## Reusable products, figures and tables

- **View manifests** (tracked): `output/views/<view>/<source>.json` for `full`, `merge`,
  `solo` and `full-pooled-label` over both sources plus `solo-deepest` over tenv only —
  nine products. They pin the recipe, members, row/condition counts and payload digests;
  the multi-GB payloads and freezes remain gitignored products in the materialised checkout.
- **Frozen truth** (view-independent): `output/truth/<source>/primary/`,
  `condition_mean/` and `diagnostics.json` — reused by every view, byte-identical on a
  re-run.
- **Ticket-03 discarded-member evidence** (tracked):
  [log/ticket-03-discarded-members.csv](log/ticket-03-discarded-members.csv) — 55,531
  per-condition rows (`solo` tenv 1,377, `solo` tahoe_bulk 52,777, `solo-deepest` tenv
  1,377) regenerated by `execution/evidence.py` from the published `members.csv` companions,
  never from a payload matrix; `execution/tests/test_evidence.py` reconciles the table with
  the released population constants.
- **Report tables** under `output/report/`:
  [metrics_primary.csv](output/report/metrics_primary.csv),
  [metrics_condition_mean.csv](output/report/metrics_condition_mean.csv),
  [metrics_per_dose.csv](output/report/metrics_per_dose.csv),
  [diagnostics_beside_scores.csv](output/report/diagnostics_beside_scores.csv),
  [paired_vs_full.csv](output/report/paired_vs_full.csv),
  [paired_vs_full_solo_deepest.csv](output/report/paired_vs_full_solo_deepest.csv),
  [baseline_mean_cell.csv](output/report/baseline_mean_cell.csv),
  [fold_spread.csv](output/report/fold_spread.csv),
  [seed_spread.csv](output/report/seed_spread.csv),
  [learning_curves.csv](output/report/learning_curves.csv),
  [run_status.csv](output/report/run_status.csv), with digests in
  [report_manifest.json](output/report/report_manifest.json).
- **CSV-only figure data and renderer**: `figure/data/*.csv` and
  `figure/script/render.py` (a thin wrapper over `execution/render.py`); the renderer needs
  no matrices, checkpoints or truth and re-renders every PNG byte-identically.
- **No weights are published.** `output/runs/` holds the accepted run evidence but no
  `output/model/`; the checkpoints stay gitignored.
- **Figures** (each rendered from the figure CSVs alone):
  - [learning_curves_tenv.png](figure/learning_curves_tenv.png) — per-view validation MSE
    by epoch, one panel per fold, both seeds.
  - [metrics_tenv.png](figure/metrics_tenv.png) — equal-fold mean ΔPCC and ΔRMSE per view,
    with fold SD and seed SD as separate error bars.
  - [metrics_per_dose_tenv.png](figure/metrics_per_dose_tenv.png) — the per-held-dose
    scores with each fold's dose and test-profile count in its title.
  - [learning_curves_tahoe_bulk.png](figure/learning_curves_tahoe_bulk.png),
    [metrics_tahoe_bulk.png](figure/metrics_tahoe_bulk.png),
    [metrics_per_dose_tahoe_bulk.png](figure/metrics_per_dose_tahoe_bulk.png) — the Tahoe
    suite.
  - [views_effect.png](figure/views_effect.png) — paired `view − full` differences with a
    fold+seed interval, `solo-deepest` in its own panel.
- **Sealed acceptance evidence**: [log/acceptance/receipt.json](log/acceptance/receipt.json)
  (schema `tsg-view-acceptance-v1`), with the health transcripts under
  `log/acceptance/health/`.

## Reproduction

Use the `arc` interpreter; the base `python3` cannot reproduce the sealed digests. Run one
heavy job at a time under the repository caps (`10 GiB` memory, `8 GiB` swap, 10-second
telemetry), and never beside a fit.

```bash
PY=/home/chang/miniconda3/envs/arc/bin/python
cd /home/chang/arc/stage2/tsg_view

# 1. Build the nine view products (four views on each source plus the tenv-only
#    solo-deepest: 4 + 4 + 1), one view/source at a time.
#    (tickets 01-04 recorded these in log/build-01.json .. log/build-04.json)
$PY execution/build.py --view full --source tenv --output output/views \
    --staging execution/checkpoints --repo-root /home/chang/arc --record log/build-01.json
#    repeat for full on tahoe_bulk; merge/solo/full-pooled-label on both sources; and
#    solo-deepest on tenv alone (there is no solo-deepest tahoe_bulk product)

# 2. Freeze each product into the Step 04 input format.
$PY execution/freeze.py --view full --source tenv --views output/views \
    --destination execution/freezes/full/tenv --repo-root /home/chang/arc
#    repeat per view/source; the freeze verifier checks the sealed partition mapping

# 3. Freeze the view-independent truth and diagnostics (both sources).
$PY execution/truth.py --source all --output output/truth --repo-root /home/chang/arc

# 4. Seal and run the 45-fit matrix, then verify it.
$PY execution/formal.py plan   --candidate execution/checkpoints/formal
$PY execution/formal.py run    --candidate execution/checkpoints/formal
$PY execution/formal.py verify --candidate execution/checkpoints/formal

# 5. Report tables and figures (two-pass publication: tenv first, then all).
$PY execution/report.py --runs output/runs --truth output/truth \
    --report output/report --figure-data figure/data
$PY figure/script/render.py --data figure/data --output figure --suite tenv
$PY figure/script/render.py --data figure/data --output figure --suite all

# 6. Tests (fixture-only tier everywhere; real-artifact tests find the source/payloads
#    or skip explicitly).
$PY -m unittest discover -s execution/tests -p test_acceptance.py -v
$PY -m unittest discover -s execution/tests -p test_docs.py -v
$PY -m unittest discover -s execution/tests -p 'test_*.py' -v

# 7. Health checks, from the repo root (a worktree lacks the gitignored payloads).
cd /home/chang/arc
bash script/check_paths.sh
bash script/check_tracked.sh
bash script/verify_all.sh          # alone: heavy, ~52 GB hashed, 9.7 GiB page cache

# 8. Rebuild the sealed receipt from the materialised checkout.
cd /home/chang/arc/stage2/tsg_view
$PY execution/acceptance.py receipt --root /home/chang/arc --phase final --verify-binary \
    --preregistration log/acceptance/pre_registration.md \
    --health-dir log/acceptance/health --output log/acceptance/receipt.json
```

The last command hashes the accepted payloads and the 49 GB Steps 01–06 artifact inventory,
so run it alone. A worktree can build and test, but the accepted payloads, the
`execution/freezes/` and `output/runs/` trees exist only in the materialised checkout.
