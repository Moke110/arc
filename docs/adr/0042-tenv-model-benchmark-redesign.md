# Rebuild the Stage 2 tenv model comparison after the corpus refactor

**Status:** proposed (2026-11-05); requirements recorded, design interview open.

The refactored dataset corpus (ADR 0033) separates native-axis source evidence
from consumer-owned model views. The next model comparison must consume that
corpus and requalify its models rather than inherit the old seven-model registry
or reuse its scores. This record captures the new request; unresolved choices
below are not defaults or authorization to start the benchmark.

## Requirements established by the request

- Preserve the existing `stage2/03_models_cmbench/` evidence and seal the old
  comparison before replacing it with the tenv comparison. The requested new
  name is `03_modes_tenv`; the archive location is resolved in round 3.
- Repeat model discovery using Biomni and primary upstream sources. Screen all
  downloaded models against the same conditions: baseline transcriptome plus
  drug input; post-treatment transcriptome or transcriptome-change output;
  approximately 18,000 measured genes at both input and output with limited
  expansion on this workstation; public model, training and testing code.
- Exclude Squidiff, PerturbNet and XPert from the new comparison. Do not rewrite
  historical eligibility or results. In particular, the prior Squidiff record
  establishes non-convergence and an opportunity-cost exclusion, not a proven
  architectural gene-space ceiling.
- Retain ADR 0030 provenance layout for admitted upstream repositories in
  `pert_model/`. Newly admitted models without setup evidence require a run card
  and an isolated environment in `stage2/01_setup/` (ADR 0021).
- Freeze three folds and three training seeds shared across models. Each model
  produces nine training-and-test run samples, one for every fold × seed pair.
- Report AUPRC_positive@100, AUPRC_negative@100, Recall_positive@P80,
  Recall_negative@P80, PCC, Spearman, RMSE and R².
- Deliver an eight-panel metric bar figure with model bars, run-sample dots and
  error bars, plus per-model positive/negative precision curves in red/blue with
  light bands showing the run-sample range. Round 3 adds both simple baselines
  to these figures and requests a third, best-single-run PR overview.

## Decisions resolved in interview round 1

- **Cold-drug evaluation.** Hold every dose and replicate of a test drug out
  of training and validation; put stereoisomers sharing the same active-parent
  connectivity into one split block. Round 3 replaces the originally agreed
  structure-derived drug input with one-hot identity; the held-out-drug split
  remains agreed, but its scientific interpretation requires clarification.
- **Complete CM gene axis.** All models consume and predict the same ordered
  18,164 genes. Input/output dimension and batch-size changes are allowed;
  HVG substitution, zero-filled reconstruction and a new core architecture are
  not admission paths.
- **Publication gate retained.** Exclude preprint-only models, retaining ADR
  0025's publication-status restriction. Public code alone does not establish
  admission.
- **Directional truth sets.** For each drug–dose group's true Δ, positive truth
  is the largest 100 strictly positive entries; negative truth is the smallest
  100 strictly negative entries. Keep fewer when fewer same-sign genes exist.
  Rank all genes by predicted Δ for positive and by its negation for negative.
  AUPRC@100 is average precision with those labels, not a truncation of the
  prediction ranking. Recall@P80 is the maximum recall attained at precision
  at least 0.8 with a nonempty prediction set, or zero if none is attainable.
- **Cache and evaluate separately.** Generate each fold's truth sets once and
  reuse them across seeds and models. Persist every `(model, fold, seed)`
  prediction before downstream evaluation. After the prediction stage, compute
  metrics in bounded CPU parallelism; reuse each direction's ranking/PR
  computation for AP, Recall@P80 and the curve rather than sorting repeatedly.
  Seeds never alter test truth.
- **PR figure.** Horizontal axis is recall and vertical axis is precision.
  Positive/negative are red/blue; the nine-run mean is the central curve and
  their minimum–maximum range is the light band. This is a PR curve, not a
  curve of AUPRC against rank.

## Decisions resolved in interview round 2

- **Fold–seed design.** Freeze three outer drug-structure-block test folds;
  each drug is tested exactly once. Within each fold's remaining approximately
  two thirds of drugs, hold out 15% for validation, preserving blocks. Approximate
  overall train/validation/test fractions are 56.7%/10%/33.3%; integer block
  membership determines exact counts. Split seed is `20261105`, training seeds
  are `17`, `42`, `2026`. Every model uses identical membership for each fold,
  independent of its training seed.
- **Consumer view and baseline.** Project tenv raw counts onto the ordered CM
  gene axis, then recompute CP10K and log1p on that axis. Each treated well's
  baseline input is the mean log-expression of same-plate DMSO wells. Test-plate
  DMSO is an allowed measured input; test treated expression is evaluator-only
  ground truth and must not be exposed to training or prediction.
- **Target and scoring unit.** Score Δ = treated log-expression minus matched
  DMSO mean. Convert post-treatment-expression outputs to that same Δ. Average
  well-level predicted and observed Δ within each `(drug_key, dose)` group,
  compute all eight metrics per group, then take an equal-weight group macro
  average for each run sample.
- **Compound and dose scope.** Use all 460 structure-resolved tenv small
  molecules, without a pharmacology-qualification filter. Exclude blanks and
  non-DMSO controls from treated training targets. All entrants must support
  actual dose as an input; an architecture's configurable gene dimension alone
  does not establish qualification.
- **Initialization and budget.** Train all formal comparisons from random
  initialization with no external pretraining. Follow author training methods,
  use validation for checkpoint selection and perform no additional HPO.
  Each `(model, fold, seed)` has a four-hour training ceiling (up to 36 hours
  per model). Mark timeouts and failures explicitly; never impute missing run
  samples. Models without nine completed runs cannot receive a complete ranking.
- **Dedicated scheduling agent.** During execution one agent owns dispatch and
  monitors WSL available RAM and GPU memory at high frequency. The user's
  instruction permits concurrent training when measured resources allow it,
  including long fits; this refines the default serial-long-fit rule for this
  experiment. Retain the repo-wide maximum of two heavy jobs, combined peak
  approximately 28 GiB, per-process memory caps, 12 GiB available-RAM admission
  floor, bounded data access and peak accounting. GPU/host profiles must justify
  starting a second model; concurrency is opportunistic, not guaranteed.
- **Figures.** For each metric the model bar is the arithmetic mean of its
  nine run samples, dots are those samples and the error bar is ± one sample
  standard deviation. This describes run variation, not independent biological
  replicates. For PR, first macro-average groups within each run on a shared
  recall grid; then compute the nine-run mean and minimum–maximum band.
- **Undefined metrics.** A direction with no true genes has NaN AP and recall;
  undefined PCC/Spearman and zero-variance-truth R² also remain NaN. Macro means
  use valid groups only and report valid-group counts for every metric.

## Decisions resolved in interview round 3

- **Model boundary.** The user defines a model by the parameters updateable in
  one backward pass. A drug encoder trained independently on data outside the
  perturbation dataset is a separate model and is not included in this test.
- **One-hot drug inputs.** Replace structure features with one-hot drug identity
  for every entrant, using its required drug-representation dimension. This
  explicitly replaces round 1's chemical-structure-input admission gate. Exact
  dimension handling where a model expects fewer than 460 coordinates is open;
  no collisions, hashing or non-one-hot encoding are implicitly allowed.
- **Author validation loss.** Select checkpoints using each author's own loss
  on the validation data, matching its training objective and reducing redundant
  metric work. Do not substitute a common Δ-MSE selection criterion. The
  across-model comparison still uses the eight common test metrics.
- **Baselines in all figures.** Include `no_change` (Δ = zero) and `delta_mean`
  (training-group mean Δ) in both the metric figure and PR figure, evaluated
  under the same folds and scoring rules as entrants. These deterministic
  baselines have no trainable model parameters; seed copies are identical and
  are not independent fits.
- **Third figure.** Add a single PR overview with one curve per model/baseline,
  each from that entry's best `(fold, seed)` run. The best-run selection score
  and how to combine directions in one curve remain open.
- **Scheduling and retries.** Sample and record RAM/GPU usage every five
  seconds. Admit at most two heavy tasks using measured peaks plus 20% headroom,
  still satisfying the approximately 28 GiB combined budget and at least
  12 GiB available RAM. Stop dispatch below the floor; if pressure continues
  worsening, save and stop the newer task belonging to this experiment.
  Repairable environmental failures permit one retry; OOM permits one retry
  after reducing batch size within the same cumulative four-hour training
  budget. Persistent failures are reported rather than replaced with synthetic
  results.
- **Archive and active paths.** Archive the old step at
  `stage2/_archive/03_models_cmbench/`, preserving its results and existing
  `_frozen` inputs. Use `stage2/03_modes_tenv/` for the new comparison.

## Remaining design frontier

Clarify one-hot dimensionality, the held-out-identity task interpretation, and
the best-run PR overview's selection/direction rules. Model discovery and source
verification remain research prerequisites; no new model is declared qualified
by this proposed record. Historical structure-based exclusion reasons alone
cannot establish exclusion from the new one-hot task.

The old ADR 0027 contract remains the contract of the frozen historical study;
its five-fold single-seed entrant list and k=200 primary metric must not silently
govern the new experiment. Final choices will be recorded here as they resolve.

## See

- [Dataset corpus contract](0033-stage2-dataset-corpus.md)
- [Historical benchmark contract](0027-stage2-benchmark-contract.md)
- [External-data and step layout](0030-external-data-provenance-and-step-layout.md)
- [Historical model publication gate](0025-stage1-stage2-directory-split.md)
- [Stage 2 overview](../../stage2/README.md)
- [Execution environment](../agents/execution-environment.md)
