# Stage 2 gene-axis — public myocardial disease transcriptome datasets

**Status:** research note (read-only survey, 2026-11-13). Not an ADR; it authorises no download.
Every accession below was verified on 2026-11-13 against NCBI E-utilities (`esearch`/`esummary`,
`db=gds`) and the GEO series-matrix headers / FTP `suppl/` listings; no large file was
downloaded. Numbers in the "n (case/control)" column come from `!Sample_characteristics_ch1`
of the series matrix, not from paper abstracts. Anything not verified this way is marked
**UNVERIFIED**.

**Question.** Which public myocardial disease transcriptome series should feed the disease
component of the narrowed Stage 2 "gene axis" (replacing the frozen 18,164-gene axis), covering
HF, MI/MIRI, chemotherapy cardiotoxicity, and any heart disease with prevalence >5%?

**Bottom line.** Chemo-cardiotoxicity is already fully covered locally (repo `pert_data/`);
MIRI is already covered locally (repo `disease_data/`, animal only — human acute-MI myocardium
bulk data does not exist). The genuinely missing, high-value, cheap addition is **human
heart-failure LV tissue** plus **one human AF series** and, optionally, **one valvular/HHD
series**. Five public series (<100 MB total) are enough; the two largest HF series are
optional replication.

## Recommendation

1. **GSE141910** — human LV RNA-seq, MAGNet: **200 HF / 166 non-failing**; open, 60 MB
   per-sample CSVs. Primary HF anchor (largest verified human HF LV transcriptome).
2. **GSE116250** — human LV RNA-seq: **50 HF (37 DCM, 13 ICM) / 14 NF**; open, **7.5 MB**
   gene-level RPKM with both ENSG and symbol columns. Cheapest HF anchor, and its **ICM arm is
   the only human post-MI-proxy LV transcriptome** we found.
3. **GSE243674** — human **iPSC-CM** anthracyclines (DOX, DNR, EPI, MTX, trastuzumab + vehicle;
   6 female donors; 3 h/24 h): **72 wells**; open, 2.3 MB counts. Best new disease-side series
   because the model (hiPSC-CM) matches the Stage 2 cell model.
4. **GSE128188** — human atrial appendage RNA-seq: **5 AF / 5 SR, paired LA+RA**; open, 0.6 MB
   counts. Covers the >5% AF category.
5. **GSE282618** — human heart RNA-seq, hypertension vs control: **6 vs 4 samples**; open,
   0.7 MB counts. Covers the >5% hypertensive-heart-disease category, but it is tiny — treat as
   exploratory only.
6. **Do not download** for MIRI or anthracycline cardiotoxicity: repo `disease_data/` already
   has GSE160516 (mouse MIRI), GSE108940 (rat MIRI), GSE137293 (pig scar-vs-border, severity
   gradient only), and `pert_data/` already has GSE217424 (DToxS iPSC-CM) and GSE244740
   (ScreenSeq hiPSC-CM cardiotox panel), both of which contain doxorubicin/daunorubicin/
   epirubicin/idarubicin/mitoxantrone. See `## Already-local assets`.

**Smallest viable first pass** if downloads are constrained:
**GSE116250 (7.5 MB) + GSE243674 (2.3 MB) + GSE128188 (0.6 MB) ≈ 10.4 MB**, combined with the
already-local MIRI + DToxS/GSE244740 anthracycline data. This covers HF (etiology-resolved,
incl. human ICM), chemo-cardiotoxicity (new independent human iPSC-CM series) and AF; HHD,
valvular AS, microarray HF replication and human acute MI stay open. Total first-pass disk
cost is 10.4 MB.

**Optional replication tier (not needed for a first pass):** GSE57338 (HF microarray, 313
samples, 1.3 GB CEL), GSE5406 (HF microarray, 210 samples, series matrix only), GSE36961 (HCM),
GSE262690 (aortic stenosis, 915 MB RDS), GSE79768 (AF microarray, 125 MB CEL). Do not fetch
these before the axis-selection rule is fixed.

## Candidate table

Legend: "open?" = downloadable with no controlled access (GEO public FTP); "ready DEG?" = a
ready-made DEG table ships with the series (else DEGs must be computed); "map?" = how cleanly
the matrix maps to human symbols. All `n` values are from the series matrix, not abstracts.

| # | Accession (GEO) | Disease | Species | Tissue / model | Human? | n (case / control) | Platform | open? | Supp files (~size) | Ready DEG? | Map to human symbols | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [GSE141910](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910) | HF (DCM 166, HCM 28, PPCM 6) | human | LV free wall, bulk | yes | **200 HF / 166 NF** | RNA-seq, GPL16791 (HiSeq 2500) | yes | `GSE141910_RAW.tar` 60.0 MB (366 per-sample ENSG CSVs) + `filelist.txt` | no — compute | ENSG IDs, one symbol-map step | **include — HF anchor** |
| 2 | [GSE116250](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116250) | HF (DCM 37, ICM 13) | human | LV, bulk | yes | **50 HF / 14 NF** | RNA-seq, GPL16791 | yes | `GSE116250_rpkm.txt.gz` **7.5 MB** | no — compute (RPKM) | `Gene` ENSG + `Common_name` symbol — cleanest of the set | **include — HF + human ICM proxy** |
| 3 | [GSE57338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57338) | HF (ischemic 95, IDCM 82) | human | LV, bulk | yes | **177 HF / 136 NF** | microarray HuGene-1.0-st, GPL11532 | yes | `GSE57338_RAW.tar` **1.3 GB** (CEL) + `filelist.txt`; series matrix also exists | no — compute (RMA/limma) | probe → symbol (GPL11532 annotation) | optional — HF replication / platform diversity |
| 4 | [GSE5406](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5406) | HF (ICM 108, IDCM 86) | human | LV, bulk | yes | **194 HF / 16 donor** | microarray HG-U133A, GPL96 | yes | none — processed **series matrix only** (size not checked) | no — compute | probe → symbol (GPL96) | optional — independent ICM cohort; only 16 controls |
| 5 | [GSE36961](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36961) | HCM (cardiomyopathy) | human | myectomy LV, bulk | yes | **106 HCM / 36 control** | microarray GPL15389 | yes | `GSE36961_non-normalized.txt.gz` 37 MB | no — compute | probe → symbol | optional — HCM >5% only if cardiomyopathy is counted separately |
| 6 | [GSE1145](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1145) | HF, mixed etiologies (IDCM, ICM, valvular, alcoholic, HCM) | human | LV, bulk | yes | 107 total, per-group n not enumerated in matrix | microarray GPL8300 + GPL570 | yes | `GSE1145_RAW.tar` **665 MB** (CEL); **no series matrix (404)** | no — compute | probe → symbol, two platforms | skip — large, no processed matrix, etiology labels coarse |
| 7 | [GSE160516](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160516) | MIRI | mouse | LV, bulk | no | 12 I/R (4 each 6/24/72 h) / 4 sham | Affymetrix Clariom S, GPL23040 | yes | already local (`source/GSE160516_deg_results.csv`) | **yes** (limma, tracked) | mouse→human ortholog (repo `map_to_human.R`) | **local — MIRI primary reference** |
| 8 | [GSE108940](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108940) | MIRI | rat | LV (above/below LAD), bulk | no | 12 I/R / 12 sham (WT + MuRF1-Tg arms) | Agilent GE2-v5 8×60K | yes | already local (`source/GSE108940_WT_deg_results.csv`) | **yes** (limma, tracked) | rat→human ortholog (repo) | **local — MIRI** |
| 9 | [GSE137293](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137293) | post-MI scar | pig | scar vs border zone, bulk | no | 8 pigs × 2 regions (16 libraries) | Illumina NovaSeq | yes | already local (`source/GSE137293_COUNTS.txt.gz`, `*_deg_results.csv`) | **partial** (Welch, tracked) | pig→human ortholog (repo) | **local, but scar-vs-border = severity gradient, not disease-vs-healthy** |
| 10 | *(no human acute-MI bulk series found)* | MI / MIRI | human | — | — | — | — | — | — | — | — | **no public transcriptome series found**; human MI exists only as single-cell/spatial: [Zenodo record 6578617](https://zenodo.org/records/6578617) ("Raw Data – Part 2: Spatial multi-omic map of human myocardial infarction"), a different evidence type |
| 11 | [GSE243674](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243674) | chemo-cardiotoxicity (TOP2 inhibitors) | human | **hiPSC-CM**, 6 healthy female donors, 0.5 µM, 3 h / 24 h | yes | 72 wells (12 vehicle; 12 per drug) | RNA-seq, GPL21697 | yes | `GSE243674_Counts_RNA_ERMatthews.txt.gz` **2.3 MB** (ENTREZID × 72) | no — compute (ENTREZ → symbol) | Entrez → symbol map needed (pitfall) | **include — best new chemo series** |
| 12 | [GSE242692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242692) | chemo-cardiotoxicity (doxorubicin) | human | donor-specific hiPSC-CM, DOX 0/0.3/1/3 µM | yes | 88 wells (donors × 4 doses; 40 "resistant"/48 "sensitive") | RNA-seq, GPL24676 | yes | `GSE242692_countsMatrix_raw.txt.gz` 2.1 MB; `GSE242692_RAW.tar` 20 MB | no — compute | gene symbols | optional — donor-variability replication |
| 13 | [GSE313746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313746) | chemo-cardiotoxicity (doxorubicin) | human | iPSC-CM, 6 donors, DOX 0.5 µM, tx0/24/144 h | yes | 42 wells (21 DOX / 21 vehicle) | RNA-seq, GPL34295 | yes | `GSE313746_Counts_RNA_EMPfortmiller.txt.gz` 1.5 MB | no — compute | gene symbols | optional — long-recovery time course |
| 14 | [GSE128188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128188) | atrial fibrillation | human | paired LA + RA appendage | yes | **5 AF / 5 SR** (20 libraries) | RNA-seq, GPL18573 | yes | `GSE128188_CountsV2_EdgeR.txt.gz` **0.6 MB** | no — compute (EdgeR counts) | gene symbols | **include — AF** |
| 15 | [GSE79768](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79768) | atrial fibrillation | human | LA + RA appendage | yes | **14 AF / 12 SR** (26 arrays) | microarray HG-U133 Plus 2.0, GPL570 | yes | `GSE79768_RAW.tar` **125 MB** (CEL) | no — compute | probe → symbol | optional — larger AF, but CEL-only |
| 16 | [GSE262690](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262690) | valvular / pressure overload (aortic stenosis) | human | LV, bulk | yes | **22 AS / 36 control** | RNA-seq, GPL24676 | yes | `GSE262690_vandy_AS_full_processed.forDE.rds.gz` **915 MB** | RDS "for DE" object (needs R) | depends on RDS object | optional — only AS-vs-control series found; big RDS |
| 17 | [GSE236191](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236191) | aortic stenosis ± T2D | human | LV myocardial biopsy | yes | **47 AS (30 non-diabetic / 17 T2D), no healthy control** | RNA-seq, GPL18573 | yes | `GSE236191_RAW.tar` 11 MB + `filelist.txt` | no | gene symbols / counts | skip — no disease-vs-healthy contrast |
| 18 | [GSE282618](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282618) | hypertensive heart disease | human | LV + RV, 5 subjects | yes | **6 hypertension / 4 control** | RNA-seq, GPL30173 | yes | `GSE282618_counts.txt.gz` 0.7 MB; VST 2.1 MB | no — compute | gene symbols | **include (exploratory) — HHD** |

Notes on the table:

- GSE141910, GSE116250, GSE57338, GSE5406, GSE36961, GSE1145 all have a GEO series matrix
  (verified by fetching its header); GSE137293 has no series matrix (404) and ships a counts
  table instead, consistent with the repo's `author_log.json`.
- "Ready DEG?" is "yes" only where the repo already tracked a limma/Welch output
  (`disease_data/*/source/*_deg_results.csv`); no *new* candidate series ships a ready DEG
  table — GEO series matrices/files are expression matrices, so the DEG step must be run
  in-repo (limma / DESeq2), which is also what the axis needs (effect sizes, not just gene
  lists).
- GSE243674's count matrix is keyed by **ENTREZID**, not symbol — budget one mapping step
  (HGNC/Entrez table already pinned under `pert_data/Annotation/`).

## Disease-prevalence taxonomy (criterion d)

There is no single authoritative "fraction of all heart disease" table, so the taxonomy below
is built from the **GBD 2019 cause hierarchy** and the headline prevalence counts published in
its JACC update ([Roth et al., *J Am Coll Cardiol* 2020;76:2982–3021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7755038/)).
The denominator is "diseases of the heart" only: GBD's **IHD + HHD + RHD + non-rheumatic
valvular disease (NRFVD) + AF/AFL + cardiomyopathy/myocarditis + congenital heart anomalies +
endocarditis**, deliberately **excluding** stroke (101 M prevalent), PAD (113 M) and aortic
aneurysm, which are vascular not myocardial disease.

Verified GBD 2019 headline values (same source): IHD **197 M** prevalent cases, 9.14 M deaths;
AF/AFL **59.7 M**; RHD **40.5 M**; degenerative mitral-valve disease **24.2 M** (a *subset* of
NRFVD — so the valve share below is a lower bound); HHD **18.6 M**; congenital heart anomalies
**13.3 M**; endocarditis 1.09 M incident cases; cardiomyopathy/myocarditis 9.14 M DALYs
(no global prevalence count is given in the paper text).

Approximate share of heart disease (arithmetic on the values above; denominator
197 + 59.7 + 40.5 + 24.2 + 18.6 + 13.3 ≈ **353 M**, before cardiomyopathy and endocarditis,
which are unquantified here):

| Category | Share of heart-disease prevalence | >5%? | Public transcriptome series mapped (≤2) |
|---|---:|---|---|
| Ischaemic heart disease / CHD (incl. MI, ICM) | ≈56% | **yes** | human LV: **GSE116250 ICM arm** (13 ICM/14 NF), **GSE5406** (108 ICM/16 donor), **GSE57338** (95 ischemic); animal MIRI already local (GSE160516, GSE108940) |
| Atrial fibrillation / flutter | ≈17% | **yes** | **GSE128188** (5 AF/5 SR, processed counts); GSE79768 (14 AF/12 SR, CEL only) |
| Rheumatic heart disease | ≈11% | **yes** (globally; <5% in US/EU) | **no public transcriptome series found** — GEO search (2026-11-13) returned only small ncRNA / valve-interstitial-cell studies (e.g. [GSE185581](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185581), [GSE88803](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88803)); no RHD-vs-healthy myocardium series |
| Non-rheumatic valvular disease (AS/MR) | ≥6.8% (MV subset only) | **yes** | **GSE262690** (22 AS/36 control, 915 MB RDS); GSE236191 exists but has **no healthy controls**; valve-tissue series (e.g. GSE153555) are not myocardium |
| Hypertensive heart disease | ≈5.3% | **yes** (borderline) | **GSE282618** (6 hypertension/4 control; very small) |
| Cardiomyopathy / myocarditis (DCM, HCM, PPCM, alcoholic, viral) | unquantified in GBD text (9.14 M DALYs) | unclear | DCM/HCM arms of GSE141910 and GSE116250; GSE36961 (106 HCM/36 control); **no adequate bulk myocarditis series** (only small ICI-myocarditis / spatial studies) |
| Congenital heart anomalies | ≈3.8% | no | none usable for adult LV myocardium |
| Endocarditis | <1% | no | no transcriptome series found |

Two caveats a reviewer will (rightly) raise:

1. **GBD prevalence ≠ incidence ≠ DALYs.** On deaths instead of prevalence, IHD is an even
   larger share and HHD/RHD smaller; the >5% set is stable for IHD, AF and NRFVD, and
   borderline for HHD/RHD. State which metric is used if this becomes an ADR.
2. **Geography matters.** RHD is >5% globally but rare in US/EU high-income populations; if
   the programme targets high-income drug markets, RHD drops out and valve/HHD/AF stay. A
   US-specific AHA breakdown would be the alternative taxonomy; it was not retrievable in
   machine-readable form in this environment (**UNVERIFIED**), so GBD is used.

## Pitfalls and caveats

- **L1000 landmark-only filtering is the biggest trap already in the repo.** Every
  `disease_data/*/local/deg/*.json` is a human-symbol **L1000 978-landmark** list (documented
  in `disease_data/README.md` and each `author_log.json`). They must **not** be reused to build
  the new axis: 978 genes is a hard ceiling and would defeat the purpose of replacing the
  18,164-gene axis. Recompute DEGs on the full transcriptome from the series matrices and
  intersect with `tenv` HVGs.
- **Species ortholog mapping** for all MIRI data (rat, mouse, pig): the repo's
  `map_to_human.R` uses homologene via Entrez and drops one-to-many mappings; one-to-many and
  missing orthologs bias the animal contribution, and pig annotation is the weakest
  (GSE137293 uses direct symbol + homologene).
- **Sample size vs DEG-count inflation.** GSE57338 (n=313) and GSE141910 (n=366) will produce
  far more DEGs than GSE128188 (n=10) or GSE282618 (n=10) at the same FDR. Rank genes by
  effect size × cross-series replication, not by DEG count; otherwise the axis is a
  large-cohort-bias axis.
- **Disease-tissue vs cell-model mismatch.** GSE141910/GSE116250 are bulk LV tissue containing
  fibroblasts, immune cells and vasculature at disease-dependent fractions; Stage 2 predicts a
  Δ in a pure hiPSC-CM monolayer. Cell-composition shifts can masquerade as CM-intrinsic
  disease signal. Either deconvolve or restrict to CM-expressed genes (e.g. intersect with
  hiPSC-CM-expressed genes from `tenv`).
- **Etiology and batch heterogeneity.** MAGNet (GSE141910) mixes DCM/HCM/PPCM; GSE116250
  separates DCM vs ICM; GSE57338/GSE5406 mix ICM vs IDCM. Pooling without etiology covariates
  or a stratified analysis makes the "HF axis" an average over mechanistically different
  diseases. Donor hearts are also not healthy hearts (brain death, ICU stay, cold cardioplegia
  — stated explicitly in the GSE141910 series summary).
- **Scar-vs-border is not disease-vs-healthy.** GSE137293 (pig) compares infarct scar with
  border zone inside the same animal — a severity gradient. It should not be interpreted as a
  disease-vs-control contrast in the axis.
- **Processed-data traps.** GSE243674 is ENTREZID-keyed; GSE262690 ships only a 915 MB
  "forDE" RDS (R required); GSE5406 has no supplementary file at all (series matrix only);
  GSE57338 is CEL-only in `suppl/`; GSE1145 has no series matrix (404) — each adds a
  different preprocessing cost, and mixing them means harmonising platforms.
- **No human acute-MI bulk data.** The negative GEO search is a real result: human MI
  transcriptomes are single-cell/spatial (e.g. the Kuppe et al. map, Zenodo record 6578617),
  not bulk case–control. Human ICM explants (GSE116250, GSE5406) are the closest human proxy
  but reflect end-stage remodelled heart, not acute ischaemia–reperfusion.
- **Disease prevalence is not the same as relevance to drug-induced Δ.** AF/AS/HHD tissues
  reflect years of remodelling; whether they belong in a *drug-response* axis is a scientific
  choice, not just a coverage choice (see open questions).

## Already-local assets

Checked so nothing is downloaded twice (paths relative to repo root):

| Local item | What it is | Verified content relevant here |
|---|---|---|
| `disease_data/GSE160516/` | mouse MIRI (I/R vs sham), limma; the workflow's Primary rank reference | 16 samples: 4 sham + 12 I/R (4 each 6/24/72 h); tracked `source/GSE160516_deg_results.csv`; `local/deg/GSE160516_MIRI.json` (L1000-filtered) |
| `disease_data/GSE108940/` | rat MIRI vs sham, limma | 24 samples (WT + MuRF1-Tg, sham/IRI, above/below LAD); tracked DEG csv; `local/deg/*.json` (L1000-filtered) |
| `disease_data/GSE137293/` | pig scar vs border zone | 16 libraries (8 pigs × scar/border); **severity gradient, not disease-vs-healthy** |
| `disease_data/*/local/deg/*.json` | 3 human-symbol DEG lists | **all L1000-978-only** — do not reuse for the new axis |
| `pert_data/GSE217424/` (DToxS-iPSC) | iPSC-CM drug panel, GSE217421 Conv + GSE217423 DGE | 2.66 GB; series-matrix `drug name` values include **Doxorubicin, Daunorubicin, Epirubicin, Idarubicin** (also Afatinib, Amiodarone, Axitinib, Bortezomib, Imatinib, Lapatinib, …) → criterion (c) already local |
| `pert_data/LINCS-LDS-1593/` | DToxS LINCS L1000, primary human CM | 872 MB; `Metadata/Small_Molecule_Metadata.txt` contains Doxorubicin, Daunorubicin, Epirubicin, Idarubicin, Mitoxantrone (verified by grep) |
| `pert_data/GSE244740/` | Rosell-Hidalgo ScreenSeq hiPSC-CM cardiotoxicity panel | 37.5 MB; `treatment` values include **Doxorubicin, Idarubicin, Mitoxantrone, 5-Fluorouracil, Bortezomib, Lapatinib, Imatinib, Dasatinib, Cyclophosphamide** → second local anthracycline source |
| `pert_data/GSE262419/` | `tenv` TempO-Seq hiPSC-CM panel | 4,206 wells, 459 drug parents × 3 doses, 18,164-gene axis (`stage2/01_bases/README.md`) |
| `pert_data/GSE247541/` | Xu COVID-19 antiviral hiPSC-CM panel | out of the Stage 2 union scope — do not add |
| `pert_data/Annotation/` | pinned HGNC complete set | symbol mapping table for Entrez-keyed series (GSE243674) |

**Consequence:** criterion (c) chemotherapy cardiotoxicity does **not** need a new download for
coverage; a new human iPSC-CM series (GSE243674) is only worth adding as an *independent,
disease-side* reference. MIRI (criterion b) is covered by three local animal series; the only
gap is human MI, which does not exist as bulk case–control data.

## Open questions for the user

1. **Axis-selection rule.** Is the rule "top-N disease DEGs by |effect| × replication,
   intersected with tenv HVG", or a fixed gene set (e.g. union of disease DEGs)? The dataset
   set is only defensible once this is fixed — e.g. GSE128188 at n=10 may contribute zero
   genes if a replication filter is used.
2. **Are AF / valvular / HHD in scope at all?** They pass the >5% prevalence test but their
   tissue reflects chronic remodelling, not a drug-perturbable state; including them may add
   chronic-disease signature rather than drug-response signal.
3. **Tissue vs cell model.** Should disease DEGs come from bulk LV tissue (GSE141910 /
   GSE116250) or only from hiPSC-CM systems (GSE243674 + local DToxS)? This determines
   whether HF tissue data enters the axis.
4. **DEG vs direction.** Stage 3 consumes a predicted Δ and a phenotype direction. Does the
   disease component need signed disease direction (up/down in disease), or only gene
   membership? GMT-style lists are insufficient for a signed axis.
5. **How much download is acceptable?** A 10.4 MB first pass is trivially safe; the optional
   tier includes 915 MB (GSE262690 RDS) and 1.3 GB (GSE57338 CEL). Confirm before fetching.
6. **Do we need raw counts or is a planned contrast on processed values acceptable?** This
   decides whether GSE57338/GSE5406 (microarray, processed matrix) are usable at all in this
   pipeline.
7. **RHD geography.** Is the programme US/EU-centric (RHD <5%, exclude) or global (RHD >5%,
   but no usable public transcriptome series anyway)?
