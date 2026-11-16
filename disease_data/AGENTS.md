# AGENTS.md — disease_data

Twenty-three external 项 (ADR 0030). Each holds `script/` (fetch/build/verify recipes + the
shared R helpers), `source/` (upstream GEO bytes, gitignored and immutable),
`local/` (declared additions: the full-transcriptome DEG CSV + tracked DEG JSON, and
the normalized log-expression pair where the series ships one), tracked
`author_log.json`, `manifest.sha256` over `source/`, and an item `.gitignore`.
Category overview and conventions: `README.md`. Layout contract: ADR 0030.

**Tiered admission to the gene axis (ADR 0071, superseding ADR 0066's
exclusion).** Whole-tissue bulk is no longer ineligible by rule: ADR 0071 admits
human atria/ventricle bulk at **tier 3**, entered only after tiers 1–2 fall short
of 2,000 axis-resolved DEGs, with the tier recorded. The published tier-3 sources
are `GSE165303`, `GSE116250` and `GSE141910` (`nHF`) and `GSE269882` (`AF`);
`GSE128188` is an owned tier-3 candidate that is not declared (it was not needed
once AF's tier-1 sources were measured). Four items remain outside the axis for
reasons that are not tissue class: `GSE282618` (human HHD bulk, 3 vs 2 exploratory
subjects — never admitted), `GSE137293` (pig infarct border zone, supplementary
evidence only) and the two mouse whole-LV MIRI series `GSE160516` / `GSE108940`
(Stage-1 L1000 rank references, not human gene-axis DEGs). Do not describe the
six whole-tissue human items as a single ineligible block; only `GSE282618` is
still unexploited, and only because its design is too small.
- `GSE108940/` — mouse MIRI (*Mus musculus*, GPL7202 Agilent 4×44K); limma WT I/R vs Sham; mouse→human homologene `inTax=10090`; L1000-978 Stage-1-only; provenance corrected in #176 (adds TP53 UP).
- `GSE116250/` — human LV HF 37 DCM + 13 ICM vs 14 NF (50/14); limma on log2(RPKM+1); full DEG + expression; **tier-3** `nHF` source under ADR 0071 (43 axis-resolved v2 DEGs); also the post-MI (ICM) proxy.
- `GSE128188/` — human AF 5 vs 5 SR (paired LA+RA, 20 libraries); DESeq2 on donor-pooled counts; full DEG + expression; owned **tier-3** AF candidate, not declared under ADR 0071 (AF stopped after its tier-1 sources).
- `GSE137293/` — pig infarct scar vs border zone (severity gradient, not disease-vs-healthy); Welch t-test; supplementary only, never an axis source.
- `GSE141910/` — human LV HF 200 (DCM/HCM/PPCM) vs 166 NF; limma on deposited voom abundance; full DEG + expression; **tier-3** `nHF` source under ADR 0071 (100 axis-resolved v2 DEGs).
- `GSE160516/` — mouse MIRI (Affymetrix Clariom S); limma pooled I/R vs Sham; Stage 1 primary rank reference; L1000-978 Stage-1-only.
- `GSE243674/` — human hiPSC-CM anthracyclines + vehicle, 72 wells; DESeq2 10 drug×time-vs-vehicle contrasts; full DEG + expression; chemo-cardiotoxicity.
- `GSE282618/` — human HHD 3 vs 2 subjects (exploratory); DESeq2 on subject-pooled counts; full DEG + expression; borderline >5% HHD category.
- `GSE165303/` — human **left-ventricle bulk** RNA-seq, 101 libraries = 50 DCM vs 51 non-failing (the series' RNA arm; the ticket prose "101 DCM" was wrong); deposited log2(normalised counts) with a fixed floor, modelled as-is with limma `lmFit`; full DEG + expression; **tier-3** `nHF` source under ADR 0071 (whole-tissue bulk enters only after tiers 1-2 fall short of 2,000 DEGs, which supersedes the ADR 0066 exclusion for tier 3).
- `GSE270788/` — human LV 10x Multiome GEX SubSeries (Amrute 2024; SuperSeries GSE218392); author CM gate `cell.type.l1 == "Cardiomyocyte"` (13,803 CM / 23 donors, 52–1,576 per donor); per-donor CM pseudobulk then one joint `~0 + arm` limma-voom fit (NICM/ICM/AMI vs the shared Donor controls); feeds `nHF` (NICM) and `ICM` (ICM+AMI). Erratum: deposit arm sizes are Donor 7 / NICM 6 / ICM 5 / AMI 5, not the ticket's 6/7 (#187 comment).
- `GSE255612/` — human left-atrium 10x snRNA (18 AF / 16 control donors); author CM gate `cell_type == "Cardiomyocytes"` (32,438 CM); the deposit is log1p(CP10K) normalised expression, so the per-donor mean-CP10K pseudobulk is tested with limma (documented linear model, not a count model); feeds `AF`.
- `GSE138252/` — human atrial **PCM1+ FANS-sorted cardiomyocyte-nuclei** bulk RNA-seq (CM-nuclear bulk, not whole tissue). The 48 GSMs map to 38 library columns: the 2020-batch matrix has 18 atrial libraries (LA 5 AF + 5 control, RA 5 AF + 3 control), GSM4103429 has no expression row, and the 2021 addendum re-sequences four 2020 hearts plus 6 HF-only libraries; the 2785 matrix is a 6-column technical-replicate subset. Published contrast = left atrium, 2020 batch only, 5 AF vs 5 control (the AF arm is AF+HF vs non-failing heart; recorded). One live-implementation limma-voom fit; **1 v2 DEG** (SLC9A7) under the live worker (the audit's TMM alternative gives 2, logged as a sensitivity); **tier-1** `AF` source under ADR 0071.
- `GSE238242/` — human left-atrial-appendage 10x **multiome snRNA arm only** (7 RNA of the 14 libraries; the ~130 MB snATAC arm is not fetched); author cell table gates `cell_type == "CM"` (3,049 of 11,986 nuclei) and the per-donor sum-of-raw-counts pseudobulk (7 donors: 4 SR / 3 persistent AF) is fitted with the one live limma-voom implementation; the honest outcome is **0 v2 DEGs** (minimum padj 0.231); **tier-1** `AF` source under ADR 0071.
- `GSE269882/` — human right-atrial-appendage **whole-tissue bulk** RNA-seq, 27 libraries; the deposited series matrix gives sinus rhythm 9 / post-operative AF 8 / paroxysmal AF 5 / persistent AF 5 (the earlier audit's 10/7 split was wrong). The count matrix's patient-number columns are joined to the titles by patient number and verified with XIST/Y sex probes (26/27 agree). Three staged AF-vs-SR contrasts in one live limma-voom fit; **0 v2 DEGs** in every contrast (pooled AF-vs-SR also 0); **tier-3** `AF` source only under ADR 0071, which supersedes the ADR 0070 whole-tissue exclusion.
- `GSE35108/` — patient-derived DCM iPSC-CM vs normal iPSC-CM, one-colour Affymetrix GPL6244 RMA log2 array. The deposit has 17 libraries but only **8 are cardiomyocytes** (a 2x2 genotype x adenoviral Serca2a design, 2 per cell); the 9 fibroblast/hESC/iPSC lines are excluded, and the probes collapse to one max-mean probe per gene (19,036 genes). Patient-derived, not mutation-modelled: **tier-2** `nHF` source under ADR 0071.
- `GSE224129/` — human MYBPC3 (cMyBP-C) HCM engineered cardiac tissue, isogenic +/+ / +/- / -/-, 4 ECT libraries each; RSEM expected counts integerised by the shared helper; `~0 + group` with the two mutant arms collapsed by maximum effect; feeds `genetic_cardiomyopathy`.
- `GSE234141/` — human TTN truncating-variant DCM hiPSC-CM vs isogenic WT/WT in 2D monolayer and EHT, atrial- and ventricular-like, 3 replicates per cell (24 libraries); one joint fit over 8 composite strata with four paired TTNtv-vs-WT contrasts collapsed by maximum effect; feeds `genetic_cardiomyopathy`.
- `GSE304575/` — human LMNA p.His222Pro DCM hiPSC-CM (6) vs CRISPR-corrected isogenic control (6); single `~0 + group` contrast; the cleanest isogenic deposit, one constant batch value; feeds `genetic_cardiomyopathy`.
- `GSE305648/` — human RBM20 DCM hiPSC cardioids on WTC11, KO / PL / RQ / wildtype (3 each) with the library batch rotated A/B/C; `~0 + group + batch` with three mutant-vs-WT contrasts collapsed by maximum effect; feeds `genetic_cardiomyopathy`.
- `GSE253223/` — human wild-type hiPSC-CM, acquired PKP2 siRNA knockdown (4/day) vs non-targeting siRNA (1/day), days 2/4/6/8; `~0 + group + day` with the days pooled and day blocked; an ARVC-loss **knockdown**, not a patient variant (recorded); feeds `genetic_cardiomyopathy`.
- `GSE194297/` — human DMD point-mutation cardiac organoids (day 56) vs the CRISPR/Cas-corrected isogenic line (3/3); single `~0 + group` contrast; feeds `genetic_cardiomyopathy`.
- `Reichart-CM/` — non-GEO cellxgene CM-only object (collection e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb, 311,418 cells / 79 donors), named by object identity; donors split by `Primary.Genetic.Diagnosis`; genotype-negative DCM donors feed `nHF`, and the tracked axis-aligned per-donor pseudobulk + genotype table is the ticket #188 `genetic_cardiomyopathy` input (single NCCM donor excluded). Ticket #188 adds `script/build_genetic_deg.py`, which derives the genotype-positive-DCM and ARVC contrasts from that tracked export read-only (no pseudobulk re-derivation), with the derived table in `local/deg/reichart-cm_genetic_deg.csv` and a tracked `.json` summary.
- `Kuppe-MI/` — non-GEO cellxgene object (collection `8191c283-0816-424b-9b61-c3e1d6258a77`, version `4780eb7d-4c7c-4630-8705-471f0de8e27a`; 191,795 nuclei / 20 donors / 28,380 features, left-ventricle snRNA), named by object identity; author CM gate `cell_type_original == "Cardiomyocyte"` (64,510 nuclei) and pseudobulk per (donor, zone) from the authors' `major_labl` (CTRL/IZ/BZ/RZ/FZ, 25 profiles); the 500-CM-nuclei profile floor leaves 18 profiles and one joint `~0 + zone` limma-voom fit with `duplicateCorrelation(block = donor)` publishes the `IZ-CTRL` contrast (4 IZ donors vs the 4 normal donors; 3,099 v2 DEGs genome-wide, 1,646 on the frozen axis under the pinned arc limma); the other zone contrasts stay recorded artifacts. **tier-1** `ICM` source under ADR 0071; the deposit carries no time-since-infarct field, so acute-vs-chronic is recorded as undeclared.

The human series that go through the in-item DEG build follow the gene-axis interface:
`local/deg/<GSE>_deg.csv` (full transcriptome, signed effect + BH padj + case/control n;
gitignored, regenerable) and the tracked `local/deg/<GSE>_deg.json` summary, plus
`local/expr/<GSE>_log_expr.tsv.gz` and `local/expr/<GSE>_samples.csv` where the series
ships one. The three Stage-1 MIRI L1000 lists (`local/deg/*.json`) are unchanged in role and
Stage-1-only; never an axis source. `GSE160516`/`GSE108940` contribute their
full-transcriptome `source/*_deg_results.csv` (DEG-only) instead.

Pinned R: human builds with `/home/chang/miniconda3/envs/biomni_e1/bin/Rscript`; the
GSE108940 mouse remap with `R_HOME=/home/chang/miniconda3/envs/arc/lib/R
/home/chang/miniconda3/envs/arc/lib/R/bin/exec/R` (the `arc` env's `R` wrapper is mis-patched).
`script/verify_all.sh` verifies every item offline. The five new series were selected by
`docs/research/stage2-stage3-gene-axis-datasets.md`; their rank-aggregation evidence is
`docs/research/stage2-stage3-gene-axis-rank-aggregation.md`.

- `GSE106297/`, `GSE157282/` — human stem-cell CM DOX/vehicle FPKM sources; ADR 0072 TOX admission, log2(FPKM+1) non-count analysis.
- `GSE232331/`, `GSE235470/`, `GSE242692/`, `GSE255505/`, `GSE313746/` — ADR 0072 TOX tier-2 count sources; source designs and controls in each `author_log.json`.
