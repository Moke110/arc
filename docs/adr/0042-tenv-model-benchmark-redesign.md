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
