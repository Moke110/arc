# Establish a Morgan baseline for Stage 3 phenotype direction prediction

**Status: accepted scope; remaining design decisions under discussion (2026-11-06)**

Stage 3 ultimately maps a drug-induced transcriptomic perturbation Δ to cardiomyocyte phenotype direction. While Stage 2 transcriptome prediction is unfinished, `stage3/02_base_morgan` will establish a Morgan-fingerprint baseline. The user accepted a shallow MLP with one nonlinear hidden layer; Morgan and future Δ predictors share the hidden-layer and output-head architecture, with the first layer accommodating their different input dimensions.

The target is decrease, within-band, or increase in the phenotype readout relative to its declared control, rather than a therapeutic-benefit label. The initial phenotype is `cell_viability`, using the magnitude directions established by ADR 0044. The output supplies three direction probabilities through Softmax. The meaning of multiple output indicators remains to be clarified before freezing head grouping.

The baseline is designed to support a later comparison of structural information with predicted transcriptomic Δ. Freeze shared sample-splitting, labeling and evaluation protocols during baseline design; a paired comparison on jointly covered samples will be needed if future Δ coverage is incomplete. This avoids making model comparisons depend on different evaluation populations. Exact splits, context/dose handling, optimization and evaluation choices are not yet settled.
