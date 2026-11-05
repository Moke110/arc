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
  name is `03_modes_tenv`; final spelling and archive location remain open.
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
  light bands showing the run-sample range.

## Decisions resolved in interview round 1

- **Cold-drug evaluation.** Hold every dose and replicate of a test drug out
  of training and validation; put stereoisomers sharing the same active-parent
  connectivity into one split block. Admission requires a structure-derived path
  for scoring an unseen drug.
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
  Seeds never alter test truth. Degenerate empty truth-set handling is still
  to be specified.
- **PR figure.** Horizontal axis is recall and vertical axis is precision.
  Positive/negative are red/blue; the nine-run mean is the central curve and
  their minimum–maximum range is the light band. This is a PR curve, not a
  curve of AUPRC against rank.

## Remaining design frontier

Open choices include validation allocation inside three folds, seed values, baseline/control
availability, normalization and aggregation, initialization and training budgets,
failure handling, uncertainty conventions, and archive paths. No new model is
declared qualified by this proposed record.

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
