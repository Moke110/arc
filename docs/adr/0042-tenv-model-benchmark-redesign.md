# Rebuild the Stage 2 tenv model comparison after the corpus refactor

**Status:** accepted (2026-11-05). The user confirmed the complete experiment
design, exact Morgan parameters and the artifact-CLI testing boundary.
Research and admission evidence remain required before formal execution.

**Spec:** [GitHub issue #40](https://github.com/Moke110/arc/issues/40), published
with ready-for-agent. The user approved the ticket breakdown; seven fixed
implementation tickets [#41](https://github.com/Moke110/arc/issues/41) through
[#47](https://github.com/Moke110/arc/issues/47) are published as native sub-issues
with verified blocked_by edges. Further per-model tickets will be instantiated
after discovery, under the approved one-model-per-ticket scope, and attached as
native blockers of formal training before dispatch. No implementation or
training has started in this planning handoff.

The dataset corpus refactor (ADR 0033) separates native-axis source evidence
from consumer-owned views. Rebuild the tenv comparison using newly screened
models and shared chemical features, preserving the historical study rather
than inheriting its entrant registry or scores.

## Boundary and model admission

- Archive the existing step at `stage2/_archive/03_models_cmbench/`, preserving
  results and its existing `_frozen` inputs. The active replacement is
  `stage2/03_modes_tenv/`, using the user's requested spelling. Do not delete
  the frozen-input bridge while other concluded studies still depend on it.
- Repeat discovery with Biomni and primary upstream sources; screen all
  downloaded models under the same final contract. Admission requires public
  model/training/testing code, formal publication, baseline transcriptome plus
  drug and dose input, post-treatment transcriptome or Δ output, and full
  18,164-gene input/output within this workstation's measured resources.
- Exclude Squidiff, PerturbNet and XPert as explicitly requested. Preserve old
  eligibility and results. Historical Squidiff evidence establishes failed
  convergence/opportunity-cost exclusion, not an architectural gene ceiling.
- Allow input/output dimension changes and the agreed chemical adapter, retaining
  core model architecture. No HVG substitution, zero-filled full-space output,
  test-treated-expression inputs or new core architecture as admission shortcuts.
  Verify that baseline expression actually participates in prediction; freely
  learned training-observation latents alone do not demonstrate that capability.
- Follow ADR 0030 for upstream repositories in `pert_model/`; leave source
  trees untouched and declare local additions. Newly admitted models require
  native setup evidence, isolated environments and run cards in `01_setup/`.
  Full-space bounded training/prediction evidence is an additional admission gate.

## Shared chemical representation and model boundary

- **Selected by the user:** one common Morgan fingerprint representation plus
  a trainable linear dimension adapter for each predictor. One-hot and KPGT
  were considered and withdrawn; neither is the current input policy.
- Use Morgan radius 2, 2,048-bit binary vectors, atom-based invariants
  (ECFP-style, not feature-based FCFP), chirality enabled. These exact parameters
  are confirmed. Fix the RDKit version and structure-normalization
  rule; cache one canonical parent-structure fingerprint per drug and share it
  across models, folds and seeds.
- Feed the cached fingerprint through `Linear(2048, d_model)` at the existing
  drug-representation interface. Initialize the adapter and perturbation
  predictor randomly for each run and train them together only on training
  responses. Keep actual dose as an explicit model input. Record adapter
  dimensions and trainable parameter counts and disclose native-recipe deviations.
- The user's model boundary is the parameters updateable in one backward pass.
  The jointly trained adapter belongs to its perturbation predictor. A drug
  encoder trained independently outside perturbation data is a separate model,
  excluded from this comparison; Morgan extraction has no learned parameters.
- Common input information does not imply identical learned adapted vectors,
  or reproduction of each author's complete original chemical pipeline.
  Models without a compatible representation interface must be reviewed rather
  than given an unapproved replacement backbone.

## Input, target and fold–seed table

- Build a consumer-owned view from `02_dataset_corpus/` tenv raw counts:
  project to the ordered 18,164-gene CM axis, then recompute CP10K and log1p.
  Retain all 460 structure-resolved small molecules without a pharmacology
  filter; blanks and non-DMSO controls are not treated training targets.
- The baseline for each treated well is its same-plate DMSO mean log-expression.
  Test-plate DMSO is an allowed input; test treated expression is evaluator-only.
  Δ is treated log-expression minus that matched-control mean. Convert any
  post-expression prediction to this common Δ.
- Average observed and predicted well-level Δ within each `(drug_key, dose)`
  group. Compute metrics per group and give groups equal macro-average weight.
- Freeze three outer cold-drug test folds, with every drug tested exactly once.
  Keep all doses/repeats and same-parent-connectivity stereoisomers together.
  Split metadata may use chemistry; test responses never determine membership.
- Within each outer fold's remaining approximately two thirds of drugs, hold
  out 15% for validation, preserving structure blocks. Overall proportions
  are approximately 56.7%/10%/33.3%; integer block sizes determine exact counts.
- Split seed is `20261105`; training seeds are `17`, `42`, `2026`.
  All models share the same frozen drug-role membership. One
  `(model, fold, seed)` is one run sample, giving nine per model rather than
  nine independent biological replicates.

## Training, scheduling and failures

- Random initialization, no external pretraining, author training methods and
  losses, no extra HPO. Select checkpoints by each author's own validation
  loss matching its training objective, rather than a common MSE or PR metric.
- Four-hour cumulative training ceiling per run, up to 36 hours per model.
  Mark incomplete runs/timeouts; never impute samples. Entries missing any of
  nine completed runs cannot receive a complete ranking.
- A dedicated scheduling agent records WSL available RAM and GPU usage every
  five seconds and owns dispatch. The user authorizes parallel model training
  when measured resources allow, including long fits. Retain at most two heavy
  tasks repo-wide, approximately 28 GiB combined peak budget, per-process caps,
  at least 12 GiB available RAM, bounded data access and peak accounting.
- Use measured task peaks plus 20% headroom to admit concurrency; otherwise
  run serially. Stop dispatch below the RAM floor; if pressure keeps worsening,
  save and stop the newer task belonging to this experiment. Do not stop
  unrelated processes.
- A repairable environmental failure permits one retry; OOM permits one retry
  with smaller batch size within the same four-hour training budget.
  Persistent failures remain reported failures.

## Evaluation and computation

- Report AUPRC_positive@100, AUPRC_negative@100, Recall_positive@P80,
  Recall_negative@P80, PCC, Spearman, RMSE and R², all on the common Δ target.
- For each group's true Δ, positive truth is its largest 100 strictly positive
  entries; negative truth is its smallest 100 strictly negative entries.
  Keep fewer if fewer same-sign genes exist. Rank all genes by predicted Δ
  for positive and by its negation for negative. AP@100 means average precision
  against those labels, not truncating prediction rankings to 100 genes.
- Recall@P80 is maximum recall at precision at least 0.8 with a nonempty
  prediction set, or zero if no such set is attainable.
- Empty-direction AP/recall, undefined PCC/Spearman and zero-variance-truth R²
  remain NaN. Macro means use valid groups and disclose valid-group counts.
- Cache truth/labels once per fold. Persist all model/fold/seed predictions
  before downstream test scoring; then evaluate in bounded CPU parallelism.
  Share rankings/PR computations between AP, Recall@P80 and curves. Author
  validation loss remains available during training for checkpoint selection.
- Include `no_change` (zero Δ) and `delta_mean` (training-group mean Δ) under the
  same folds and scoring rules, in all three figures. Their identical
  seed copies are deterministic references, not independent trained fits.

## User-authorized partial evaluation amendment (2026-11-06)

The user explicitly amended the evaluation timing contract: “Amend the contract;
score finished models now”. Completed runs with validated persisted predictions
may now be scored against evaluator-only test truth while the remaining formal
runs continue. This exception supersedes the global prediction-barrier requirement
for partial evaluation only; the canonical final evaluation still requires all
63 model/fold/seed outcomes to reach terminal status, including recorded failures.

- Freeze each partial evaluation's included run identities, prediction hashes,
  completion counts and truth/metric identities. Use a separate, explicitly dated
  partial-evaluation namespace for metrics, tables, source arrays and snapshot
  figures. Publish the current user-facing PNG figures at
  `stage2/03_modes_tenv/output/metrics.png`,
  `stage2/03_modes_tenv/output/pr_panels.png` and
  `stage2/03_modes_tenv/output/best_run.png`, with SVG companions when
  produced. Canonical final figures overwrite these same display paths after all
  63 outcomes settle. Interim provenance, scores and snapshots stay separate;
  display-file replacement does not replace canonical evaluation evidence.
  Partial work must neither create a final seal nor overwrite canonical final
  evaluation records or the live scheduler's prediction/status records.
- Keep training and evaluation independent: partial test results must not influence
  training, checkpoint selection, hyperparameters, retries, admission or scheduling.
  Existing author validation-loss checkpoint selection and the frozen run matrix
  remain binding. Evaluation is bounded CPU work subject to current measured RAM
  and the repo-wide heavy-job limit; ongoing fits continue under their scheduler.
- Primary partial figures compare entries with all nine completed runs and include
  both deterministic reference baselines. Show any partly completed entry separately
  with its completed/expected run count; it cannot receive a complete ranking or a
  nine-run mean/range label. Retain NaN reasons and valid-group counts. Best-single-run
  displays remain explicitly post-hoc, with their available-run selection scope.
- Partial delivery supplements the study. After all 63 outcomes settle, perform the
  canonical integrity/seal, full scoring, replay, independent review and tracker
  closure under the original final-delivery requirements. Partial evidence does
  not satisfy final matrix completeness or close formal execution by itself.

## Figures

1. Eight metric panels: model/baseline bars are nine-run arithmetic means,
   dots are individual run samples, error bars are ± one sample standard
   deviation. Describe run variation without treating runs as independent
   biological replicates.
2. Per-model/baseline PR panels: x = recall, y = precision, positive red and
   negative blue. Macro-average groups on a shared recall grid within each
   run, then plot the nine-run mean and light minimum–maximum range.
3. One best-single-run PR overview: choose each entry's run maximizing
   `(AP_positive@100 + AP_negative@100) / 2`, break ties in fixed fold/seed
   order, and average its two directional precision curves on a common recall
   grid into one line. Label this as post-hoc best test performance; entries
   may select different folds. Nine-run summaries remain the primary evidence.

## Research evidence and remaining work

Biomni Layer 2 was invoked in the prescribed environment: sandbox DNS failed;
network-enabled broad and narrow PubMed calls returned HTTP 500; an arXiv call
was cancelled without a result. A further bounded Biomni Google-search call
failed on sandbox DNS and returned HTTP 429 on its network-enabled retry;
its process has ended. These attempts are not a completed Biomni
research pass; an effective retrieval is still required.

Primary pages identify new candidates
[PertDiff](https://academic.oup.com/bioinformatics/article/42/7/btag482/8722292),
[PerturbDiff](https://proceedings.mlr.press/v306/yuan26d.html) and
[scDFM](https://openreview.net/pdf?id=QSGanMEcUV). Publication is not full-space,
adapter or resource admission evidence. Formal publication of theislab
CellFlow has not been verified. Reassess every candidate under the selected
Morgan-plus-adapter contract; no new entrant is declared qualified here.

The [RDKit Morgan documentation](https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints)
defines configurable circular fingerprints. It does not establish superiority
over KPGT on tenv. KPGT was inspected but not selected; its local provenance
gap and 2,304-dimensional base output do not block the chosen Morgan design.

The user confirmed the shared understanding and authorized the published spec
and dependency-ticket workflow. Test at the artifact CLI boundary: frozen
inputs, isolated-environment fit/predict, persisted predictions, offline metrics
and all three figures; add focused behavior checks for metric edge cases and
resource scheduling. The approved graph starts at research #41 and input #42;
#42 → #43 → #44, {#41,#44} → #45, then per-model admission tickets → #46 → #47.
ADR 0027 remains the contract of the frozen historical benchmark.

## See

- [Dataset corpus contract](0033-stage2-dataset-corpus.md)
- [Historical benchmark contract](0027-stage2-benchmark-contract.md)
- [External-data and step layout](0030-external-data-provenance-and-step-layout.md)
- [Publication gate](0025-stage1-stage2-directory-split.md)
- [Stage 2 overview](../../stage2/README.md)
- [Execution environment](../agents/execution-environment.md)

## Amendment (2026-11-06): joint CMonge and PerturbDiff v2 protocols

The user explicitly confirmed archival, joint formal migration and resume after
checks, and requested more frequent checkpoint saving, validation and early-stop
detection for **both CMonge and PerturbDiff**. This supersedes the earlier
PerturbDiff-only amendment: replace all **18 fold/seed identities** together,
nine per model; the other **45 identities and artifacts remain byte-identical**.
These are adapted validation-stopping recipes, not claims to reproduce the
complete author recipes. Keep each model's architecture, native objective,
optimizer, training population/batch, schedule and sampling settings unchanged.

### CMonge: `cmonge-v2-full-validation-fixed-draws`

- Validate and save checkpoints every **50 iterations**, plus the configured
  final iteration; retain the native **10,000-iteration** maximum. Each event
  evaluates **all 138 sorted validation conditions** (46 drugs) using identical
  seeded draws across events and an **equal-weight mean of native condition
  objectives**. Preserve training RNG state around validation. Legacy v1 sampled
  one held-out condition per event every 100 iterations; these populations differ.
- Initialize the patience reference at the first eligible event at/after **500
  iterations**. Stop after **10 subsequent** eligible checks without relative
  improvement of at least **0.001**, using denominator
  `max(abs(reference), 1e-12)`; earliest patience stop is **1,000 iterations**.
  Select the strict minimum checkpoint independently from that patience threshold,
  including earlier events. Persist event history and latest/best prediction-only
  checkpoints; these snapshots do not establish exact training continuation.
- The bounded full-axis train/save/reload/predict proof completed in **82.29 s**,
  peak RSS **1,634,500 KiB (~1.56 GiB)**, with finite 2 × 18,164 predictions and
  zero reload error. A separate 51-iteration cost probe timed out: its first
  50 training iterations took 39.0069 s and full validation took 71.5771 s.
  Linear extrapolation to 10,000 iterations is **~6.14 hours**, above the original
  four-hour allowance; earliest patience eligibility projects ~2,211 s. These
  synthetic probes prove execution/cost only, not formal qualification,
  convergence or completion. Keep the full population and report honest timeouts.

### PerturbDiff: `perturbdiff-v2-all-validation-ema`

- Validation and checkpoint eligibility occur every **1,000 steps**, plus the
  configured maximum **200,000 steps**. Deterministic EMA validation uses fixed
  draws and **all 138 validation groups**, in bounded batches of eight; the final
  two-group batch receives its group-count weight. Record IDs, count and hash.
  Legacy v1 used the first 40 groups, covering only 14 of 46 validation drugs.
- Initialize the reference and zero patience counter at the first eligible event
  at/after **20,000 steps**. Stop after **10 subsequent** eligible checks without
  an absolute reference improvement of at least **0.01**; earliest patience stop
  is **30,000 steps**. Independently select the strict minimum EMA validation-loss
  checkpoint across all eligible events, including those before 20,000.
  Legacy-history variation is operational context, not a biological noise
  estimate or confidence interval for the all-group validation population.
- The 1,001-step CUDA synthetic smoke is **TEST ONLY**, with all 138 groups,
  prediction reload and timeout rejection; it proves neither convergence nor a
  completed formal v2 run.

### Joint archival, budget and activation contract

Restart both v2 fits from their declared seeds. Exact resume is unsupported;
legacy snapshots and prediction reloads do not prove training continuation.
Preserve the original identities, native binaries/checkpoints, logs, histories
and hashes of **four completed legacy v1 outcomes** (CMonge and PerturbDiff fold
0 seeds 17/42) and **one interrupted legacy PerturbDiff fold 0 seed 2026**.
Archive them as legacy, never relabel them v2 or credit them to the 18 formal
v2 slots. Final matrix size remains 63 outcomes.

The original **14,400-second cumulative allowance per model/fold/seed** spans
v1, v2 and every retry. Conservatively debit all legacy scheduler attempt wall
time, including setup/inference when inseparable, and reconcile interrupted work
against logs/checkpoints. Preserve an authenticated, identity-bound ledger;
seed restart grants no new credit. CMonge fold 0 seed 17 has conservative prior
wall debit **13,975.509756658023 s**, leaving **424.490243341977 s**; seed 42 debit
is **9,056.932746226995 s**, leaving **5,343.067253773005 s**. Independently
reconcile every attempt and untouched slot before migration. The seed-17
remaining allowance cannot reach the earliest patience stop under measured cost.
New training includes setup, validation and checkpoint overhead; separately timed
new inference has its own watchdog. Both models require an independent hard
wall-deadline supervisor that terminates the owned process group; cooperative
between-group checks alone cannot enforce the cap. Exhaustion yields an honest
terminal timeout, never successful adoption of a partial fit.

**Execution status after independent review:** joint migration is **APPLIED**.
The protected atomic exchange committed after actual-host stopped-state guards
and independent integration, ledger/deadline and renderer-barrier review passed.
The reviewed stage inventory SHA256 is
`73fb21036a2d04a12c5d8d9ce712ed544b1c38eb70fbfff884dc6753c8b60028`;
`output/formal/joint-v2-apply.json` has SHA256
`09a5d90779a8e2dfafee58b6f2274b090ac3bb9b3d1ff9e307b632156cd6d547`.
Post-apply verification derived all 63 identities: 18 changed and 45 unchanged,
with those 45 queue records unchanged. All 22 archive operations and complete
old/new formal inventories matched their journal hashes. The 18 ledger rows
retain the original 14,400-second caps without a reset. The five formal JSON
records are tracked; the retained `execution/.formal-v2-generation/` recovery
generation and native legacy binaries remain local and untouched.

At the 2026-11-06 23:29:54 actual-host check, the owned scheduler service
was **ACTIVE** (MainPID 1487683, `Restart=no`) in startup validation. Host RAM
available was ~40 GiB and GPU usage was 1,431 MiB. Runtime hashing/validation
was in progress; first-worker dispatch and its budget account/telemetry had not
yet been verified. The committed journal's `restart_authorized: false` records
that APPLY itself did not restart training; the separately owned launch follows
the user's explicit authorization. The earlier approval rejection preceded
that confirmation; the joint APPLY passed automatic review. Neither APPLIED nor
service activation establishes a completed v2 fit.

The initial startup subsequently failed closed before worker dispatch because
`nvidia-smi` was absent from the service PATH. Its first budget account had no
attempts and retained the historical debit unchanged. The launch owner verified
the old service inactive/MainPID 0 and lock free, then restarted the same owned
unit with `/usr/lib/wsl/lib` on PATH (invocation
`9cc05786bff54bb3a037840a21316c30`). Integrity checks were running again;
first-worker/account/telemetry proof remained pending at this documentation
snapshot. No protocol, formal identity or cumulative cap changed for that fix.

Final delivery after all 63 outcomes settle must regenerate learning curves
alongside metrics/PR figures and after every replay. Both models' formal v2 and
legacy v1 histories have separate populations/panels, source hashes and native
coordinates; missing, interrupted and timeout histories stay explicit. Never
splice seed-restarted trajectories or count legacy histories as formal v2 slots.
The standalone renderer implements both v2 schemas (commit `d03f44a`); this is
presentation support, not production activation or convergence evidence.
Issues **#40, #46 and #47 remain open** pending formal outcomes, final evaluation,
replay, delivery and review. This amendment is not closure.

## Amendment (2026-11-07) — CMonge v2 fold 0 seed 17 supplemental attempt

The user authorized one supplemental attempt for **CMonge v2, fold 0, seed 17**
after the original scheduler exits cleanly. This attempt starts from seed
initialization with the frozen recipe, validation population, inputs and
scientific run identity. Its separate authorization allows **14,400 seconds
for training, including watchdog setup, and 1,800 seconds for inference**.
The original cumulative-budget timeout, its ledger and its artifacts remain
separate evidence. The original attempt had **424.49 seconds** available and
exhausted the original **14,400-second** cumulative allowance.

Review and pin the one-shot policy against the settled original run records,
then dispatch under the shared scheduler lock and the dedicated scope's
resource and wall-time limits. A durable reservation spends this authorization.
Before final scoring, independently review the terminal evidence and select
this attempt for the existing slot. Completion contributes one scored run.
An authenticated failure or timeout occupies the slot without a score.
Missing, invalid or mismatched evidence blocks final scoring. Selection is
fixed by the authorization before test metrics are computed.

Final evaluation uses a separate materialized overlay with the selected
attempt and the other 62 original outcomes. Keep both attempts' native
histories separate from one another and from legacy v1. Disclose the unequal
budget, both statuses, measured phase timings and host RSS in the final
tables and figures and on each replay.

**Executed 2026-11-07.** The original scheduler settled 62 completions and
the target timeout and exited cleanly. Independent reviewers authenticated the
policy, terminal outcome, 63-slot selection and evaluation before publication.
The supplemental run completed by validation patience and supplies exactly one
of the 63 scored slots. The original timeout, original budget and separate
history remain in the evidence. The committed [final prediction manifest](../../stage2/03_modes_tenv/execution/final_record/predictions.json)
and [evaluation provenance](../../stage2/03_modes_tenv/execution/final_record/evaluation-provenance.json)
bind that decision to the final scores. The [step results and replay command](../../stage2/03_modes_tenv/README.md#final-results)
record the resource measurements, historical scheduler limits, defined-run
counts and final figure locations.
