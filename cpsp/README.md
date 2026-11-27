# CarPS-Phenotype (CPSP)

CarPS-Phenotype aims to predict the direction of change in one cardiomyocyte
phenotype under a specified condition. The completed work establishes the
Heliste cell-viability corpora and a structure-and-dose comparator for a later
predictor using CarPS-Transcriptome responses.

## 01_corpus — frozen Heliste viability measurements

Heliste 2021 Dataset S1 reports ATP-dependent CellTiter-Glo luminescence in
HL-1 mouse atrial cardiomyocytes at 0.04, 0.2, 1 and 5 µM. The release keeps
duration-matched normoxia (`CTRL`, baseline) and hypoxia–reoxygenation
(`HYP`, H-R) as separate `cell_viability` corpora. H-R used 24 hours of pretreatment, 48 hours
of hypoxia and 24 hours of reoxygenation. Both arms ran for about 96 hours.
Each arm retains 2,744 source measurements spanning 686 labels mapped to 663
desalted parents, with 12 measurements pending identity resolution. The release
preserves the original readouts and source aliases and leaves replicate
identity unspecified. [Corpus report](01_corpus/README.md) ·
[Corpus manifest](01_corpus/output/manifest.json) ·
[Phenotype source archive](../pheno_data/README.md)

The derived magnitude classes are −1 below 80, 0 from 80 through 120 inclusive
and +1 above 120. They use the working assumption that the same-arm 0.01% DMSO
vehicle reference is 100. The exact oxygen-specific normalization denominator
remains unresolved. [Corpus decision](../docs/adr/0044-stage3-source-first-corpora.md)

## 02_base_morgan — completed structure-and-dose benchmark

Morgan uses a chirality-aware 2,048-bit fingerprint of the desalted parent
and log dose standardized on the training partition. The formal benchmark
completed all 60 fits across the two arms, three training seeds and ten fixed
partitions (three cold-drug, four cold-dose and three warm-start folds).
Cold drug had the lowest mean macro-F1 in both arms, at 0.4065 for baseline and
0.4115 for H-R, and every cold-drug fit had zero recall for the +1 class. The
independent audit reproduced the metrics and aggregates from 49,392 saved test
predictions. The migrated release replayed all 60 fits from saved artifacts.
[Formal results](02_base_morgan/README.md#formal-results-and-release-evidence) ·
[Release record](02_base_morgan/execution/release_record.json) ·
[Independent audit](02_base_morgan/execution/validation/formal-data-audit.md)

## Prospective CarPS-Transcriptome interface

The planned predictor takes signed CarPS-Transcriptome Δ on the selected
ordered gene axis, concatenated with dose. Each condition–phenotype pair has
its own model with the Morgan baseline's 128-unit ReLU hidden layer and
three-class softmax, returning probabilities in the order `[-1, 0, +1]`.
The paired evaluation will reuse the frozen measurement identities, labels,
partition roles and metrics on rows covered by both predictors. The CPST
prediction source for the Heliste comparison, coverage and integration remain
pending. [CarPS-Transcriptome configuration](../cpst/README.md#selected-carps-transcriptome-configuration) ·
[Comparison protocol](../docs/adr/0045-stage3-morgan-direction-baseline.md) ·
[Frozen splits](02_base_morgan/execution/split_manifest.json)

For build, verification and replay instructions, see the [CPSP guide](AGENTS.md).
