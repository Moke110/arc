# Denoise-DEG prediction axes

Denoise-DEG selects genes for cardiomyocyte perturbation prediction by
screening for expression detection and DMSO variability, then ranking retained
genes by their presence in tenv perturbation and disease DEG lists. The accepted
release publishes five nested subsets of the fixed 18,164-gene CM axis, from
`denoise-1k` to `denoise-5k`. [Gene ranking](output/release/gene_ranking.csv) ·
[Release manifest](output/release/manifest.json)

## From screening to ordered axes

1. **Screen expression and DMSO variability.** The native tenv corpus provided
   4,478 nonblank wells, including 192 DMSO wells across 16 plates. A gene
   failed the detection screen if its count was nonzero in fewer than 10% of
   nonblank wells. For the variability screen, a one-way random plate
   intercept fitted by moments separated between-plate SD from within-plate
   residual SD on well-level `log2(1 + CPM)` using native library depths. For
   each SD, the boundary slope was the 90th percentile of SD divided by
   mean DMSO expression on that scale among the 18,159 DMSO-detected genes.
   Values strictly above either boundary qualified for exclusion. The three
   rules together excluded 3,093 genes and retained 15,071.
   [Denoising panels](figure/denoise.png) ·
   [Per-gene exclusions](figure/data/diagnostics_gene_exclusions.csv) ·
   [Population and corpus pin](figure/data/diagnostics_population.json)
2. **Count tenv and disease DEG presence.** Plate-wise limma-voom processing
   compared measured parent–dose conditions with matched DMSO. Within each
   contrast, hits required Benjamini–Hochberg-adjusted `p < 0.05` and
   `|log2FC| > 1.5`. The highest `|log2FC| × −log10(padj)` resolved repeated
   observations of a gene for the same parent–dose across plates. Up to 100
   hits per direction and parent–dose contributed to `p_tenv`, yielding
   77,348 votes. Each of the four eligible diseases contributed its published
   tier-first DEG1000 to `p_disease`. [Score summary](figure/data/diagnostics_score_summary.json) ·
   [Step 02.1 disease evidence](../02.1_hvg-deg/README.md)
3. **Rank retained genes and publish memberships.** The score is
   `p_tenv × 2^p_disease`. Of the retained genes, 13,166 had positive scores.
   Each axis takes the highest-scoring genes for its size and breaks ties
   by ascending frozen gene index. The axis CSVs store `gene_index` and
   `symbol` in ascending frozen gene-index order. [Score distributions](figure/deg-scores.png) ·
   [Axis membership by score](figure/score-in-axis.png) ·
   [Membership plot data](figure/data/score-in-axis.csv)

## Accepted evidence and downstream use

Disease coverage divides the intersection with each disease's DEG1000 by the
complete list size, including genes excluded during screening. Coverage for
Denoise-DEG 5k ranged from 43.4% for genetic cardiomyopathy to 48.2% for TOX.
The release acceptance checked all 20 disease-by-axis values, verified
nesting and projected eight native tenv rows per axis. The native DEG replay
reproduced the tenv vote counts. [Coverage plot](figure/disease-DEG1000-coverage.png) ·
[Coverage values](figure/data/disease-DEG1000-coverage.csv) ·
[Release acceptance](log/acceptance.json)

Step 03 published plate-pooled tenv training views for Denoise-DEG 2k through
5k. Step 04 compared the 3k, 4k and 5k axes with TranSiGen and selected
Denoise-DEG 5k for the subsequent studies. Step 09 compared the 2k through
5k axes with UniMol and selected `denoise-5k` for CarPS-Transcriptome.
[Step 03 training views](../../03_data_train/README.md#tenv-training-axes) ·
[Step 04 results](../../04_tsg_axis/README.md#formal-results) ·
[Step 09 results](../../09_cpst_axis/README.md#formal-results) ·
[Selection decision](../../../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md)

Use the shared `prediction_axis.open_axis("denoise-5k")` reader to verify
and load the ordered genes. The reproduction guide covers both release rebuilds
with the pinned inputs and figure replay from published plot data. The
rendering receipt records the current PNG, plot data and renderer digests.
[Prediction-axis API](../02.1_hvg-deg/output/api/prediction_axis.py) ·
[Reproduction](AGENTS.md#reproduction) ·
[Current rendering receipt](../log/342-figure-render.json)
