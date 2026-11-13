# AGENTS.md — disease_data

Eleven external 项 (ADR 0030). Each holds `script/` (fetch/build/verify recipes + the
shared R helpers), `source/` (upstream GEO bytes, gitignored and immutable),
`local/` (declared additions: the full-transcriptome DEG CSV + tracked DEG JSON, and
the normalized log-expression pair where the series ships one), tracked
`author_log.json`, `manifest.sha256` over `source/`, and an item `.gitignore`.
Category overview and conventions: `README.md`. Layout contract: ADR 0030.
- `GSE108940/` — mouse MIRI (*Mus musculus*, GPL7202 Agilent 4×44K); limma WT I/R vs Sham; mouse→human homologene `inTax=10090`; L1000-978 Stage-1-only; provenance corrected in #176 (adds TP53 UP).
- `GSE116250/` — human LV HF 37 DCM + 13 ICM vs 14 NF (50/14); limma on log2(RPKM+1); full DEG + expression; HF cheapest anchor and post-MI (ICM) proxy.
- `GSE128188/` — human AF 5 vs 5 SR (paired LA+RA, 20 libraries); DESeq2 on donor-pooled counts; full DEG + expression; the >5% AF category.
- `GSE137293/` — pig infarct scar vs border zone (severity gradient, not disease-vs-healthy); Welch t-test; supplementary only, never an axis source.
- `GSE141910/` — human LV HF 200 (DCM/HCM/PPCM) vs 166 NF; limma on deposited voom abundance; full DEG + expression; HF primary anchor.
- `GSE160516/` — mouse MIRI (Affymetrix Clariom S); limma pooled I/R vs Sham; Stage 1 primary rank reference; L1000-978 Stage-1-only.
- `GSE243674/` — human hiPSC-CM anthracyclines + vehicle, 72 wells; DESeq2 10 drug×time-vs-vehicle contrasts; full DEG + expression; chemo-cardiotoxicity.
- `GSE282618/` — human HHD 3 vs 2 subjects (exploratory); DESeq2 on subject-pooled counts; full DEG + expression; borderline >5% HHD category.
- `GSE270788/` — human LV 10x Multiome GEX SubSeries (Amrute 2024; SuperSeries GSE218392); author CM gate `cell.type.l1 == "Cardiomyocyte"` (13,803 CM / 23 donors, 52–1,576 per donor); per-donor CM pseudobulk then one joint `~0 + arm` limma-voom fit (NICM/ICM/AMI vs the shared Donor controls); feeds `nHF` (NICM) and `ICM` (ICM+AMI). Erratum: deposit arm sizes are Donor 7 / NICM 6 / ICM 5 / AMI 5, not the ticket's 6/7 (#187 comment).
- `GSE255612/` — human left-atrium 10x snRNA (18 AF / 16 control donors); author CM gate `cell_type == "Cardiomyocytes"` (32,438 CM); the deposit is log1p(CP10K) normalised expression, so the per-donor mean-CP10K pseudobulk is tested with limma (documented linear model, not a count model); feeds `AF`.
- `Reichart-CM/` — non-GEO cellxgene CM-only object (collection e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb, 311,418 cells / 79 donors), named by object identity; donors split by `Primary.Genetic.Diagnosis`; genotype-negative DCM donors feed `nHF`, and the tracked axis-aligned per-donor pseudobulk + genotype table is the ticket #188 `genetic_cardiomyopathy` input (single NCCM donor excluded).

The five human series follow the gene-axis interface: `local/deg/<GSE>_deg.csv` (full
transcriptome, signed effect + BH padj + case/control n; gitignored, regenerable) and the
tracked `local/deg/<GSE>_deg.json` summary, plus `local/expr/<GSE>_log_expr.tsv.gz` and
`local/expr/<GSE>_samples.csv`. The three Stage-1 MIRI L1000 lists (`local/deg/*.json`) are
unchanged in role and Stage-1-only; never an axis source. `GSE160516`/`GSE108940` contribute
their full-transcriptome `source/*_deg_results.csv` (DEG-only) instead.

Pinned R: human builds with `/home/chang/miniconda3/envs/biomni_e1/bin/Rscript`; the
GSE108940 mouse remap with `R_HOME=/home/chang/miniconda3/envs/arc/lib/R
/home/chang/miniconda3/envs/arc/lib/R/bin/exec/R` (the `arc` env's `R` wrapper is mis-patched).
`script/verify_all.sh` verifies every item offline. The five new series were selected by
`docs/research/stage2-stage3-gene-axis-datasets.md`; their rank-aggregation evidence is
`docs/research/stage2-stage3-gene-axis-rank-aggregation.md`.
