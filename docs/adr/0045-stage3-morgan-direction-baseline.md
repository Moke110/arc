# Establish a Morgan baseline for Stage 3 phenotype direction prediction

**Status: agreed design; awaiting final shared-understanding confirmation (2026-11-06)**

Stage 3 ultimately maps a drug-induced transcriptomic perturbation Δ to cardiomyocyte phenotype direction. While Stage 2 transcriptome prediction is unfinished, `stage3/02_base_morgan` will establish a Morgan-fingerprint baseline. The user accepted a shallow MLP with one nonlinear hidden layer; Morgan and future Δ predictors share the hidden-layer and output-head architecture, with the first layer accommodating their different input dimensions.

The target is decrease, within-band, or increase in the phenotype readout relative to its declared control, rather than a therapeutic-benefit label. The initial phenotype is `cell_viability`, using the magnitude directions established by ADR 0044. Train a separate model for each `(condition, phenotype)` pair. Each model outputs only the probabilities of -1, 0 and +1 through one three-class Softmax; these probabilities sum to one. There is no shared multi-phenotype output model. The Morgan baseline takes both the fingerprint and dose as inputs; the accepted dose encoding is specified below.

The baseline is designed to support a later comparison of structural information with predicted transcriptomic Δ. Freeze shared sample-splitting, labeling and evaluation protocols during baseline design; a paired comparison on jointly covered samples will be needed if future Δ coverage is incomplete. This avoids making model comparisons depend on different evaluation populations. Evaluate three generalization settings: unseen drugs, unseen doses, and warm start. The accepted split, encoding, training and reporting rules are specified below. These three settings expand the initial cold-drug-only recommendation.

## Generalization settings and metrics (2026-11-06)

The user accepted three distinct evaluation settings. Cold drug groups all source aliases and doses of each desalted parent together, so test parents never occur in training; dose levels may overlap. Cold dose uses four leave-one-dose-out evaluations over 0.04, 0.2, 1 and 5 µM: the held dose is absent globally from training, while test parents have observations at other training doses. Report each held dose separately and distinguish interpolation at interior doses from extrapolation at boundary doses. Warm start holds out parent–dose combinations while ensuring both the parent and dose individually occur in training; all source records for a parent–dose combination stay together.

Use macro-F1 as the primary numerical metric in every setting, with per-class precision/recall, confusion matrices and log loss for probability quality. Validation construction, repetition counts and aggregation follow the configuration below. These settings are distinct experiments, rather than three properties required of a single test split.

## Baseline configuration and repetition (2026-11-06)

Use the Stage 2 ADR 0042 training seeds `17`, `42`, `2026` for every fold, with fold membership fixed independently of training seed and shared across baseline/H-R models and future Δ comparisons. Cold drug has three outer test folds; within each outer training pool, hold out 20% of parents for validation. Cold dose has four held-dose folds, with validation parent–dose combinations drawn from the remaining doses while validation parents remain represented in training. Warm start uses three outer parent–dose test folds, with each combination tested exactly once. Within each non-test pool, hold out 20% of combinations for validation. Preserve training coverage of validation/test parents and doses, and keep all source records for each combination together. This replaces the earlier 70%/15%/15% proposal. Test partitions are reserved for final evaluation.

The accepted Morgan configuration is radius 2, 2048 binary bits, with chirality enabled. Concatenate the fingerprint with `log10(dose_um)`, standardized using training-partition statistics only. Use `Linear(input,128) → ReLU → Linear(128,3)` with a three-class Softmax probability output. Future Δ predictors concatenate Δ and dose and use the same hidden and output architecture, changing the first-layer input dimension. Do not search depth or width in 02.

Train with ordinary unweighted cross-entropy and AdamW, without oversampling, and use validation macro-F1 for early stopping. This preserves a simple baseline without intentionally changing class prevalence. Use learning rate `1e-3`, weight decay `1e-4`, batch size `64`, at most 300 epochs and early-stopping patience of 30 epochs without validation macro-F1 improvement. Restore the best checkpoint, preferring the earlier epoch for ties. No additional hyperparameter search is included.

Use split seed `20261105` and persist explicit split manifests. Each condition–phenotype model has (3 cold-drug + 4 cold-dose + 3 warm-start) folds × 3 training seeds = 30 fits; the current two conditions require 60 fits. Keep fitting sequential for long training runs under the repository execution-environment rules.

Persist predictions and metrics for every fold × seed. Report metric means and standard deviations, with each cold-dose level reported separately before its overall average. For cold drug, additionally concatenate all three held-out fold predictions within each training seed and report the resulting overall metrics. Seed/fold summaries are descriptive; training seeds are not independent biological replicates.

## Delivery boundary

The next spec and implementation cover the complete `02_base_morgan` training and evaluation pipeline. Freeze and retain the shared comparison protocol now. Concrete Δ integration is a later design step after Stage 2 determines the prediction source, gene axis and coverage; do not invent Δ vectors or make their availability a prerequisite for completing 02. Proceed to specification only after the user confirms the complete design understanding.
