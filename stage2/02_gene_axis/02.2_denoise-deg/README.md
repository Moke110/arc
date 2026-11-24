# Denoise-DEG prediction axes

Five nested axes select 1,000 through 5,000 genes from the fixed 18,164-gene Stage 2 common space. The release includes a [ranked table](output/release/gene_ranking.csv), five two-column gene-axis CSVs and a [checksum manifest](output/release/manifest.json). Axis rows are in ascending frozen gene-index order. Membership takes the highest scores first, breaking ties by ascending frozen gene index.

## Selection

- Denoising evaluates 4,478 nonblank tenv wells and 192 DMSO wells on 16 plates. Exclude genes with nonzero counts in fewer than 10% of nonblank wells. Separately estimate between-plate and within-plate residual SD on well-level `log2(1 + CPM)`, using native library depths and a one-way random plate intercept estimated by moments. For each SD, set the slope of `SD = k × mean DMSO log2(1 + CPM)` to the 90th percentile of SD divided by mean expression among the 18,159 DMSO-detected genes. Exclude values strictly above either line. The union excludes 3,093 genes and retains 15,071. [Denoising panels](figure/denoise.png) · [Per-gene exclusions](figure/data/diagnostics_gene_exclusions.csv)
- Tenv presence counts parent–dose conditions whose top-DEG contains the gene. Fit limma-voom on each plate with matched DMSO as the baseline. Within each contrast, apply BH correction and require `padj < 0.05` and `|log2FC| > 1.5`. Rank by `|log2FC| × −log10(padj)`, keep the strongest result per gene and parent–dose across repeated plates, and take at most 100 per direction. Disease presence counts eligible [02.1](../02.1_hvg-deg/README.md) categories whose tier-first top 200 contains the gene. Compute both counts from these reference lists, then apply the exclusion union before ranking with `score = p_tenv × 2^p_disease`. [Score distributions](figure/deg-scores.png)
- Coverage is the intersection of an axis with the published top 200 for each eligible disease divided by 200, including any DEG200 members removed during denoising in the denominator. [Coverage plot](figure/disease-DEG200-coverage.png) · [Exact coverage values](figure/data/disease-DEG200-coverage.csv)
- Each row of the axis membership plot uses the score-distribution bins. Pale blue gives the retained-gene count, and dark blue runs down from the bar top by the count selected in that row's axis. Gene counts use a logarithmic axis. [Axis membership plot](figure/score-in-axis.png) · [Bin counts](figure/data/score-in-axis.csv)

## Reproduction

Run from the repository root with the native tenv corpus materialized. To rebuild the selection evidence, set `ARC_R_HOME` to the R installation containing limma and use a staging directory containing the two execution scripts and the R worker. Apply the repository resource limits. The fit uses all measured doses and limma-voom. The GSE262419 publication reported the highest dose with DESeq2. [Tsai et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11691792/)

For a replay from the published evidence, copy the step to a candidate directory and run the last two commands. The renderer reads plot data from CSV.

```bash
python "$STEP/execution/diagnose.py" --repo "$PWD" --step "$STEP"
python "$STEP/execution/diagnose_scores.py" --repo "$PWD" --step "$STEP"
python "$STEP/execution/build_axes.py" --repo "$PWD" --step "$STEP"
python "$STEP/figure/script/render.py" --repo "$PWD" --step "$STEP"
```
