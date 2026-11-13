# Stage 2 gene axis — cardiomyocyte-derived Control-vs-Disease sources

**Status:** research note (read-only survey, 2026-11-13). Not an ADR; it authorises no download.
Accessions were resolved live against NCBI E-utilities (`esearch`/`esummary`/`eLink`, `db=gds`),
GEO series-matrix headers, FTP `suppl/` listings and HTTP `HEAD`; cellxgene objects at header level
only (`obs` schema, shape, `raw/X` dtype). The only payload fetched was the **1.27 MB** cell table
`GSE270788_metadata.csv.gz`. `n` per arm comes from `!Sample_characteristics_ch1`, GSM
`!Sample_title`, or an h5ad `obs` cross-tab — never from an abstract; anything unverified this way
is marked **UNVERIFIED**. Inputs: the three working reports under `tmp/gene-axis-survey/`.

**Question.** Which public transcriptomes can supply **cardiomyocyte (CM)-derived**
Control-vs-Disease evidence for the narrowed Stage 2 gene axis, under the eligibility rule: a
source qualifies only as **(a)** an in-vitro CM / organoid / engineered-tissue **bulk** DEG
contrast, or **(b)** a **CM-gated pseudobulk** from single-cell/single-nucleus data — with
priority **human scRNA > hiPSC/hESC-CM > animal scRNA**. Whole-myocardium bulk is out of scope
even when it is human and disease-relevant (see *Eligibility consequences*).

**Bottom line.** The axis keeps CM-derived evidence for almost every category: non-ischaemic HF
(NICM 7 vs 6), ICM (5 vs the same 6), acute MI (5 vs 6, author-declared), AF (18 vs 16), six
genetic cardiomyopathies (isogenic hiPSC-CM/EHT/cardioid bulk), and both cardiotoxicity classes
(GSE243674 + the already-local DToxS/ScreenSeq panels). **GSE270788 (Amrute et al., *Nature*
2024) is the one dataset supplying three categories at once** — NICM-HF, ICM and AMI — on six
shared donors, with an author CM label and fetchable per-GSM counts. HHD and PPCM have **no
eligible source at all**; ARVC patient missense, HCM-MYH7, LQT/CPVT and a sham-controlled animal
chronic-IHF series with ≥3 replicates are gaps. The rule **rejects the whole-tissue bulk anchors**
(GSE141910, GSE116250, GSE128188, GSE282618, GSE160516, GSE108940) — none is CM-derived — moving
the MIRI anchor to animal snRNA/scRNA and leaving AF on one series.

### GSE270788 — the pivot series

- **23 GEX donors**; the cell metadata `GSE270788_metadata.csv.gz` (**1.27 MB, 46,384 cells**)
  carries `condition` (Donor/HF) and `HF.etiology`: **Donor 6 / NICM 7 / ICM 5 / AMI 5**.
- CM gate: `cell.type.l1 == "Cardiomyocyte"` → **13,803 CM cells**, per-donor **52–1,576**;
  **every donor ≥52 CM**, so all 23 donors can enter a pseudobulk.
- Counts are per-GSM `*_filtered_feature_bc_matrix.h5` (**23 files, ~0.93 GB** total). The series
  `GSE270788_RAW.tar` is **80.7 GB and must not be used**.
- The three disease arms share the **same 6 Donor controls**; the AMI arm carries **no per-sample
  time-since-infarct column**, so "acute" is **author-declared**.
- It is a **multiome series (GEX + ATAC)**; **GSE218392 is its SuperSeries**. Take the GEX
  SubSeries only — do not plan ATAC payloads into the fetch.

## Coverage by disease category

`n` is donors (or wells) per arm for the primary accession in the row; alternates are named in the
same cell. **COVERED** = eligible source with ≥3/arm and a usable matrix; **PARTIAL** = eligible
but confounded, n<3, or matrix-limited; **GAP** = no eligible source. This table records the raw
survey categories; ADR 0066 merges the *acute MI* row into `ICM`, so the axis has six categories,
not seven.

| Category | Candidate accession(s) | n (case / control) | CM-label source | Matrix format + size | Verdict |
|---|---|---|---|---|---|
| nHF (non-ischaemic HF) | [GSE270788](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270788) NICM arm; [Reichart DCM h5ad](https://cellxgene.cziscience.com/collections/e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb); [GSE183852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183852) | 7 / 6; 52 / 18; 18 / 27 | `cell.type.l1`; `cell_type = cardiac muscle cell` (CM-only object); Seurat `Robj` | per-GSM h5 23 files **≈0.93 GB**; `.h5ad` 3.19 GB; counts csv 465 MB | **COVERED** |
| ICM (chronic ischaemic) | GSE270788 ICM arm; [GSE298023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298023); [GSE145154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145154) | 5 / 6; 9 / 4; 15 / 5 | author label (GSE270788 only); labels only inside 4.0 GB `.rds` **UNVERIFIED**; none → re-cluster | ≈0.93 GB; `RAW.tar` 1.8 GB + rds 4.0 GB; `RAW.tar` 1.1 GB | **COVERED**, one label-clean series |
| MI — acute | GSE270788 AMI arm | 5 / 6 | author label; **no time-since-infarct column** | ≈0.93 GB | **PARTIAL** (author-declared "acute") |
| AF | [GSE255612](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255612) | 18 / 16 | author `cell_type == "Cardiomyocytes"` (32,438 LA CM cells) | `mtx.gz` counts **957 MB** + metadata 8.4 MB | **COVERED** (atrial CM) |
| genetic_cardiomyopathy — HCM | [GSE224129](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224129); [GSE140965](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140965); [GSE262876](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262876); [GSE274124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274124); [GSE251993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251993) | 4/4 (+4); 5/arm; 4/4; 3/3; 3/3 | isogenic or non-carrier hiPSC-CM lines; genotype in series matrix | 807 K / 1.1 M / 621 K / 1.5 M / 1.6 M | **COVERED** (MYBPC3, TPM1) |
| — DCM-TTN | [GSE234141](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234141) | 3 per genotype × cell type × format (24 libs) | isogenic WT/WT sibling | counts **1.1 MB** | **COVERED** |
| — DCM-LMNA | [GSE304575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304575); [GSE101514](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101514); [GSE136252](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136252) | 6 / 6; 12 libs, per-arm **UNVERIFIED**; CM subset **UNVERIFIED** | CRISPR-corrected isogenic | 798 K; 12 MB; 548 K | **COVERED** (GSE304575 is the clean one) |
| — DCM-RBM20 | [GSE305648](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305648); [GSE329393](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329393) | 3/3/3/3; 3 vs 5 (×2 families) | isogenic WTC11 / CRISPR-corrected | 748 K + 3 DESeq2 tables; 578 K | **COVERED** |
| — ARVC | [GSE289628](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289628); [GSE248187](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248187); [GSE253223](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253223) | 4/4; 3/3; 16/4 | isogenic DSP; DSP EHT; siPKP2 knockdown | 910 K; DEG table only 1.5 M; 1.5 M | **PARTIAL** (no patient-missense bulk) |
| — DMD | [GSE194297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194297); [GSE169190](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169190); [GSE310327](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310327) | 3/3; 3/3/3; 5 iso / 17 | isogenic corrected | 335 K; 659 K; xlsx 46 M | **COVERED** (low axis priority) |
| cardiotox_TOP2 | [GSE243674](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243674); [GSE242692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242692); [GSE313746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313746); local GSE217424 + GSE244740 | 12/drug + 12 vehicle; 88 wells; 21/21; local | hiPSC-CM, treatment in series matrix | 2.3 MB (ENTREZID); 2.1 MB; 1.5 MB; 2.66 GB + 37.5 MB local | **COVERED** |
| cardiotox_other | GSE243674 trastuzumab arm; DToxS (lapatinib, imatinib, afatinib, bortezomib); GSE244740 (5-FU, dasatinib, cyclophosphamide) | as above | as above | as above | **COVERED** |
| HHD | none eligible | — | — | — | **GAP** |
| PPCM | none eligible | — | — | — | **GAP** |

Note: the requested category list has no acute-MI row; it is added because GSE270788's AMI arm is
the only human CM-derived acute-MI arm found, and the third category that series supplies.

## Domain 1 — human heart scRNA / snRNA (CM pseudobulk)

| Accession / object | Disease role | Species + model | n (case / control) | Platform | CM-label source | Deposit format + size | Licence | Verdict | Evidence |
|---|---|---|---|---|---|---|---|---|---|
| [GSE270788](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270788) (Amrute) | nHF (NICM) + ICM + acute MI | human LV, multiome | 7 NICM, 5 ICM, 5 AMI / 6 Donor (23 GEX donors) | 10x Multiome (GEX+ATAC) | author `cell.type.l1 == "Cardiomyocyte"` in shipped metadata | per-GSM `*_filtered_feature_bc_matrix.h5`, 23 files **≈0.93 GB** (`RAW.tar` 80.7 GB — do not use) | GEO public; no licence field **UNVERIFIED** | **INGEST — three arms, one control set** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270788) · [PMID 39443792](https://pubmed.ncbi.nlm.nih.gov/39443792/) |
| [GSE255612](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255612) | AF | human left atrium, snRNA | 18 AF / 16 control | 10x 3′ v3 | author `cell_type == "Cardiomyocytes"` (32,438 cells) | `mtx.gz` counts 957 MB + metadata 8.4 MB | GEO public | **INGEST** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255612) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255612/suppl/) |
| [cellxgene Reichart CMP — Cardiomyocytes](https://cellxgene.cziscience.com/collections/e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb) | nHF (DCM) | human LV/RV/apex/IVS, sc+snRNA | 52 DCM / 18 normal (+8 ARVC, 1 NCCM) | 10x 3′ v2/v3 | `cell_type == "cardiac muscle cell"`; CM-only object, 311,418 cells | `.h5ad` **3.19 GB**, CSR + `raw` counts | CC BY 4.0 (portal-level **UNVERIFIED**) | **INGEST** — exclude ARVC/NCCM; `assay` confounded with cohort | [cellxgene](https://cellxgene.cziscience.com/collections/e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb) |
| [cellxgene — myocarditis CM object](https://cellxgene.cziscience.com/collections/328d71f0-0ed7-4518-966f-be6bd0797324) | myocarditis | human LV/apex, sc+snRNA | 24 myocarditis / 18 normal | 10x 3′ v2/v3 | `author_cell_type` vCM1–vCM6 (90,273 CMs) | `.h5ad` **1.38 GB** | CC BY 4.0 (portal-level **UNVERIFIED**) | **INGEST** — shares controls with row above | [cellxgene](https://cellxgene.cziscience.com/collections/328d71f0-0ed7-4518-966f-be6bd0797324) |
| [GSE183852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183852) (Koenig) | nHF (DCM) | human LV apex, 7 scRNA + 38 snRNA | 18 DCM / 27 non-diseased | 10x | Seurat `Robj` `cell_type` — no standalone table | integrated counts csv **465 MB**; Robj 2.1 / 8.1 / 12 GB | GEO public | **INGEST (B)** — CM come from snRNA libraries | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183852) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/) |
| [GSE298023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298023) | ICM | human LV incl. scar, sn/scRNA | 9 ICM (5 + 4 scar) / 4 donor (+5 DCM) | 10x | labels only inside a 4.0 GB `.rds` — **UNVERIFIED** | per-sample 10x mtx in `RAW.tar` **1.8 GB**; rds 4.0 GB | GEO public | **INGEST (B)** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298023) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298023/suppl/) |
| [GSE145154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145154) (Wang) | ICM + DCM | human LV/RV, scRNA | 15 ICM / 5 normal (DCM 10, RM 6) | 10x | none in GEO → re-cluster | 10x mtx in `RAW.tar` **1.1 GB** | GEO public | **INGEST (C)** — no labels, region-heterogeneous | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145154) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145154/suppl/) |
| [cellxgene Kuppe All-snRNA](https://cellxgene.cziscience.com/collections/8191c283-0816-424b-9b61-c3e1d6258a77) | MI | human LV, snRNA | 16 MI / 4 normal | 10x 3′ v3 | `cell_type_original == "Cardiomyocyte"` (64,510 nuclei) | `.h5ad` **1.86 GB** | CC BY 4.0 (portal) | **CONDITIONAL** — acute-vs-chronic needs an external clinical join | [cellxgene](https://cellxgene.cziscience.com/collections/8191c283-0816-424b-9b61-c3e1d6258a77) |
| [GSE181764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181764) | HCM | human IVS myectomy, scRNA | 7 HCM / 2 normal | 10x | none → re-cluster | `RAW.tar` **287 MB** | GEO public | **REPLICATION ONLY** — 2 controls | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181764) |
| [GSE138262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138262), [GSE174691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174691) | HCM | human LV / IVS, scRNA | 16 / 0 and 9 / 0 | — | — | — | GEO public | **REJECT** — zero controls | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138262) |
| [CAREBANK / PERIHEART](https://doi.org/10.1016/j.xcrm.2024.101556) | IHF / AF (atrial) | human right atrium, snRNA | ≤4 diseased donors per arm | 10x 3′ v3 | `cardiac muscle cell` present | `.h5ad` 1.61 / 2.79 GB | cellxgene | **REJECT** — atrial CM, n too small | [doi](https://doi.org/10.1016/j.xcrm.2024.101556) |
| [GSE131778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131778), [GSE203160](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203160), [GSE121893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121893) | wrong leads | — | — | — | — | — | GEO public | **REJECT** — coronary-artery SMC (not Litviňuková); bulk not single-cell; **no CMs** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131778) |

## Domain 2 — hiPSC/hESC-CM disease bulk (and the cardiotox panels the axis needs)

The bulk-model survey deliberately excluded drug-dose-only designs; the three cardiotoxicity rows
are added here because the axis needs those categories and they are the same model system.

| Accession | Disease role / mechanism | Species + model (isogenic?) | n per arm | Platform | CM-label source | Deposit format + size | Licence | Verdict | Evidence |
|---|---|---|---|---|---|---|---|---|---|
| [GSE224129](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224129) | HCM — MYBPC3 | human hiPSC ECT; MYBPC3 +/-, -/- vs +/+ (isogenic) | 4 / 4 / 4 | bulk polyA, NovaSeq 6000 | n/a — the differentiated culture *is* CM; genotype in matrix | RSEM gene matrix **807 K** | GEO public | **INGEST** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224129) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224129/suppl/) |
| [GSE140965](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140965) | HCM — MYBPC3 LoF panel | hiPSC-CM; 4 LoF clones vs CRISPR-neg (isogenic) | 5 per arm (6 arms, 30) | HiSeq 4000 | as above | STAR gene counts **1.1 M** | GEO public | **INGEST** — batch confounded with control set A/B | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140965) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140965/suppl/) |
| [GSE251993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251993) | HCM + DCM — TPM1 | hiPSC EHT; TPM1 E62Q (HCM) and E54K (DCM) vs WT (isogenic) | 3 / 3 / 3 | HiSeq 2500 | as above | raw counts **1.6 M** | GEO public | **INGEST** — split into two contrasts | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251993) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE251nnn/GSE251993/suppl/) |
| [GSE234141](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234141) | DCM — TTNtv | hiPSC-CM atrial/ventricular; TTN p.Ser21134Trpfs*25 vs WT/WT, 2D + EHT | 3 per genotype × cell type × format (24) | DNBSEQ-G50 | as above | counts **1.1 M** | GEO public | **INGEST** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234141) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234141/suppl/) |
| [GSE329393](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329393) | DCM / LVNC — RBM20 | hiPSC-CM; 6RBM vs 6crispr (LVNC), 9RBM vs 9crispr (DCM) | 3 vs 5 per family (16) | NovaSeq 6000 | as above | raw counts **578 K** + VST 2.1 M | GEO public | **INGEST** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329393) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329393/suppl/) |
| [GSE305648](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305648) | DCM — RBM20 | hiPSC cardioids; RBM20 KO / PL / RQ vs WT, WTC11 | 3 / 3 / 3 / 3 | NextSeq 2000 | as above | counts **748 K** + 3 DESeq2 tables | GEO public | **INGEST — best organoid candidate** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305648) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305648/suppl/) |
| [GSE304575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304575) | DCM — LMNA | hiPSC-CM; LMNA p.His222Pro het vs CRISPR-corrected | 6 / 6 | NovaSeq 6000 | as above | raw counts **798 K** | GEO public | **INGEST — cleanest isogenic** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304575) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304575/suppl/) |
| [GSE262876](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262876) | HCM — MYBPC3 | hiPSC 2D CM + organoids; D389V vs non-carrier | 4 / 4 (2D); 4 / 4 (3D) | NextSeq 2000 | as above | raw counts **621 K** | GEO public | **INGEST** — healthy-donor-only; do not pool 2D/3D | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262876) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262876/suppl/) |
| [GSE274124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274124) | HCM — MYBPC3 biallelic | hiPSC-CM; c.506-1G>A + c.442G>A vs WTC11 | 3 / 3 | NovaSeq 6000 | as above | **TPM only** 1.5 M | GEO public | **INGEST** — no count model | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274124) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274124/suppl/) |
| [GSE221945](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221945) | Brugada — TBX5 | hiPSC-CM; TBX5-G145R vs isogenic + population controls + KO | 4 / 3 / 4 / 3 | HiSeq 3000 | as above | counts **1.2 M** + 6 DESeq2 tables | GEO public | **INGEST** — strongest channelopathy | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221945) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221945/suppl/) |
| [GSE289628](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289628) | ARVC — DSP | hiPSC-CM; DSP p.Arg1113X/WT vs WT (isogenic) | 4 / 4 | NextSeq 2000 | as above | featureCounts **910 K** | GEO public | **INGEST (confounded)** — AAV6 vs siRNA vehicle | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289628) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289628/suppl/) |
| [GSE248187](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248187) | ARVC — DSP | hiPSC EHT; DSP-/- vs DSP+/+ (isogenic; exclude 24 rat NRVM) | 3 / 3 human | NovaSeq 6000 | as above | **DEG table only** 1.5 M | GEO public | **INGEST (DEG-only)** — no re-analysis | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248187) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248187/suppl/) |
| [GSE253223](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253223) | ARVC — acquired PKP2 loss | hiPSC-CM WT; siPKP2 vs siNeg, d2/4/6/8 | 16 / 4 | NovaSeq 6000 | as above | counts **1.5 M** | GEO public | **INGEST** — knockdown, not a patient variant | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253223) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253223/suppl/) |
| [GSE194297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194297) | DMD | hiPSC cardiac organoids d56; DMD point mutation vs corrected | 3 / 3 | HiSeq 4000 | as above | absolute counts **335 K** | GEO public | **INGEST** (low priority) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194297) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194297/suppl/) |
| [GSE169190](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169190) | DMD | hiPSC-CM; exon-44 deletion vs reframed / exon-skipped | 3 / 3 / 3 | NextSeq 500 | as above | counts **659 K** + RPKM + 3 edgeR | GEO public | **INGEST** (low priority) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169190) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169190/suppl/) |
| [GSE310327](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310327) | DMD | hiPSC-CM; DMD19 vs isogenic DMD19-iso | 5 iso / 17 derived libs | NovaSeq X Plus | as above | raw counts xlsx **46 M** + TPM 13 M + `RAW.tar` 132 M | GEO public | **INGEST (large xlsx)** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310327) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE310nnn/GSE310327/suppl/) |
| [GSE243674](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243674) | cardiotox TOP2 + other | hiPSC-CM, 6 female donors; DOX / DNR / EPI / MTX + trastuzumab vs vehicle, 3 h/24 h | 12/drug + 12 vehicle (72 wells) | HiSeq 2500 | as above | counts **2.3 MB** (ENTREZID-keyed) | GEO public | **INGEST** — one series, both classes | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243674) |
| [GSE242692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242692) | cardiotox TOP2 | hiPSC-CM; DOX 0 / 0.3 / 1 / 3 µM | 88 wells | NovaSeq 6000 | as above | raw counts **2.1 MB** | GEO public | **INGEST (optional)** — donor variability | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242692) |
| [GSE313746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313746) | cardiotox TOP2 | hiPSC-CM; DOX 0.5 µM, tx0/24/144 h | 21 / 21 | NovaSeq 6000 | as above | counts **1.5 MB** | GEO public | **INGEST (optional)** — long recovery | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313746) |
| [GSE217424](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217424) / [GSE244740](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244740) | cardiotox both | hiPSC-CM panels, **already local** under `pert_data/` | panel wells | series matrix / ScreenSeq | as above | 2.66 GB (DToxS) + 37.5 MB (ScreenSeq) | GEO public | **LOCAL — do not re-fetch** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217424) |
| [GSE305865](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305865), [GSE240775](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240775), [GSE331321](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE331321) | HCM (MYPN / TNNI3 / non-coding) | hiPSC-CM | 2/2; 1/1/1; 2 | — | — | — | GEO public | **REJECT** — n≤2 per arm | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305865) |
| [GSE176209](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176209), [GSE210783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210783), [GSE263980](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263980) | ARVC / DCM-RBM20 / HCM | hiPSC-CM | 5/5; —; 3/3 | — | — | no matrix: `data_row_count = 0` / `RAW.tar` only / header-only | GEO public | **REJECT** — no deposited expression matrix | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176209) |
| [GSE301005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301005), [GSE271775](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271775) | ARVC / HCM | hiPSC **epicardial cells** / **cardiac fibroblasts** | 17; 6/6 | — | — | — | GEO public | **REJECT** — not cardiomyocytes | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301005) |

## Domain 3 — animal MI / IHF scRNA-snRNA

| Accession | Category | Species + model | n per timepoint (case / control) | Assay + platform | CM-label source | Deposit format + size | Licence | Verdict | Evidence |
|---|---|---|---|---|---|---|---|---|---|
| [GSE146285](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146285) | acute (1/3/14 d) | mouse C57BL/6 male, LAD I/R + separate permanent-MI 3 d | sham 1 d 4, sham 14 d 8, I/R 1 d 9, I/R 3 d 7, I/R 14 d 3, MI 3 d 5 (36) | scRNA, CEL-seq2, whole LV | none shipped → marker re-derive (`Nppa/Nppb/Myh6/Ttn`) | 36 × UMI count tsv in `RAW.tar` **17.5 MB** | GEO public | **INGEST — acute anchor** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146285) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146285/suppl/) · [PMID 33398012](https://pubmed.ncbi.nlm.nih.gov/33398012/) |
| [GSE316860](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316860) | acute (24 h) | mouse WT male, LAD 45 min + 24 h reperfusion | sham 3 / I/R + vehicle 3 (+ enpatoran 3) | snRNA, 10x | snRNA retains CM nuclei; marker call, no author table | 9 × 10x triplets in `RAW.tar` **349 MB** | GEO public | **INGEST — cleanest acute I/R** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316860) · [PMID 42447539](https://pubmed.ncbi.nlm.nih.gov/42447539/) |
| [GSE308914](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308914) | acute D1/D4/D7 **+ chronic D28** | mouse C57BL/6J male+female, permanent LAD | 6 per timepoint (D0 6 = **naive**, 3M+3F) | scRNA, 10x 3′ | author CM I/II clusters in paper; GEO ships triplets → re-derive | 30 × triplets in `RAW.tar` **5.22 GB** | GEO public | **INGEST — the only ≥3-replicate chronic arm** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308914) · [PMID 41410034](https://pubmed.ncbi.nlm.nih.gov/41410034/) |
| [GSE253767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253767) | chronic (28 d) | rat Sprague–Dawley female, permanent LAD | sham 2 / MI 2 (+ 2 hydrogel arms) | snRNA, 10x | no annotation → marker-based | 8 × UMI csv in `RAW.tar` **174 MB** | GEO public | **INGEST — n=2 only, replication** | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253767) · [suppl](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253767/suppl/) |
| [GSE185289](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185289) | chronic (28 d) | pig, LAD at P28, harvest P56 | CTL 4 samples / 3 animals vs MIP28 4 samples / 2 animals | snRNA, 10x | **author CM1–CM6** in `Cardiomyocyte_Seurat.Robj.gz` | `RAW.tar` **2.42 GB** + CM Robj **2.5 GB** | GEO public | **CONDITIONAL** — 2 infarcted animals; HF phenotype **UNVERIFIED**; naive controls | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185289) · [PMID 35666813](https://pubmed.ncbi.nlm.nih.gov/35666813/) |
| [GSE227189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227189) | acute + chronic (1/7/30 dpi) | mouse LAD vs sham | **n=1 per condition** | scRNA | annotated CELLxGENE h5ad, 194,315 cells | h5ad on CELLxGENE | cellxgene | **REJECT** — zero replicates; keep the h5ad as a CM-marker reference | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227189) |
| [GSE302892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302892), [GSE193290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193290), [GSE163956](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163956), [GSE299884](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299884) | acute + chronic | mouse/rat MI or I/R | **n=1 per timepoint/arm**; several have no sham | — | — | — | GEO public | **REJECT** — single biological replicate | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302892) |
| [GSE262430](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262430), [GSE269054](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269054), [GSE176092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176092) | acute + chronic | rat I/R ± ECM; mouse MI | — | — | — | — | GEO public | **REJECT** — no sham arm at all, or spatial-only (Visium, no suspension) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262430) |
| [GSE185265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185265), [GSE280373](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280373), [GSE306718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306718) | chronic | mouse MI / LAD | — | — | — | — | GEO public | **REJECT** — sorted non-CM populations (Tcf21 non-myocyte, CD45+, mural) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185265) |

## Eligibility consequences

The rule disqualifies whole-tissue bulk that the earlier dataset survey
([`stage2-stage3-gene-axis-datasets.md`](stage2-stage3-gene-axis-datasets.md)) had used or
recommended — **none of it is cardiomyocyte-derived**:

| Disqualified series | What it is | Why it fails |
|---|---|---|
| [GSE141910](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910) | human LV bulk RNA-seq, 200 HF / 166 NF | whole-tissue LV: CM-intrinsic change is mixed with cell-composition shift |
| [GSE116250](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116250) | human LV bulk, 50 HF / 14 NF (13 ICM) | same; its ICM arm was the only human bulk post-MI proxy |
| [GSE128188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128188) | human atrial appendage bulk, 5 AF / 5 SR | whole-tissue atrium, not CM |
| [GSE282618](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282618) | human LV+RV bulk, 3 HHD / 2 control subjects | whole-tissue, and n=5 subjects |
| [GSE203160](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203160) | human LV bulk, 8 IHD / 7 control | whole-tissue; surfaced as a bulk IHF fallback, fails the same test |
| [GSE160516](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160516) | mouse LV bulk MIRI, 4 sham / 12 I/R | whole-LV bulk — the current local MIRI anchor |
| [GSE108940](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108940) | mouse LV bulk MuRF1 I/R, 12 sham / 12 IRI | whole-LV bulk (species provenance also unresolved locally) |

Consequences: (i) the **MIRI anchor moves from bulk to animal snRNA/scRNA** (GSE146285, GSE316860,
GSE308914); (ii) **AF rests on GSE255612 alone**, because GSE128188 was its only previous series;
(iii) **HHD loses its only series** and becomes a gap; (iv) the human HF anchor becomes CM-gated
(GSE270788 + the two cellxgene CM objects) rather than MAGNet. GSE137293 (pig, scar vs border) was
never eligible either — a severity gradient inside one animal, not disease-vs-control.

## Shared-control / non-independence traps

1. **GSE270788 — three arms, one control set.** NICM (7), ICM (5) and AMI (5) are all contrasted
   against the **same 6 Donor** hearts (`HF.etiology` = Donor 6 / NICM 7 / ICM 5 / AMI 5, 23 GEX
   donors). The three contrasts are not independent; the Donor signal is counted three times.
2. **Reichart DCM object ↔ myocarditis object.** Both cellxgene collections reuse the same
   H-series healthy donors (H2–H49): the 18 "normal" in each object are partly the same people.
3. **DToxS cardiotoxicity contrasts.** Vehicle wells are shared across drug classes, so the
   doxorubicin / daunorubicin / epirubicin / idarubicin / trastuzumab contrasts are correlated
   through one control set on one plate layout.
4. *(same trap, animal side)* **GSE308914** uses one Day-0 baseline (n=6, naive rather than sham)
   for both its acute (D1/D4/D7) and chronic (D28) arms.

**Handling rule.** Declare the overlap in the dataset manifest and **never pool overlapping arms
as independent replicates** — no summed `n`, no treating shared-control contrasts as mutual
cross-validation. Model donor (or well) as the unit; where arms share controls, fit one design
(`~ arm + donor`) instead of running separate case-control tests, and report the effective n per
independent control set rather than the naive total.

## Gaps

- **HHD** — no eligible source. The only human HHD series (GSE282618, 3 vs 2 subjects) is
  whole-LV bulk; no human HHD single-cell set with a CM arm and controls was found in GEO or
  cellxgene.
- **PPCM** — zero human single-cell hits; bulk transcriptomes only.
- **ARVC patient missense** — DSP is covered (GSE289628 confounded by AAV6-vs-siRNA vehicle;
  GSE248187 DEG-table only); DSG2 (GSE176209) deposited **no expression data**
  (`data_row_count = 0`); PKP2 exists only as an acute siRNA knockdown (GSE253223), and
  GSE301005 is epicardial cells, not CM.
- **HCM-MYH7** — the only human MYH7 datasets are single-cell (GSE276210 / GSE220933) or mouse
  (GSE201755); no human bulk hiPSC-CM MYH7 series with a control arm. MYBPC3 and TPM1 are
  covered; TNNI3 is n=1/arm (GSE240775, rejected).
- **LQT/CPVT** — no human bulk hiPSC-CM disease-vs-control series; only SQT1 (GSE185687,
  healthy-donor) and Brugada (GSE221945, isogenic) are usable channelopathy arms.
- **Animal chronic IHF** — no series with (a) CMs, (b) a **sham-operated** control and (c) ≥3
  replicates/arm: GSE308914 D28 has n=3/sex but a naive control, GSE253767 is n=2, GSE185289 is
  2 pigs with an unverified HF phenotype; the rest is n=1 or bulk (GSE236374 is the bulk companion).

## Verification and provenance

**How.** Accessions were resolved on **2026-11-13** with `esearch`/`esummary`/`eLink` (`db=gds`)
for series identity and platform; arm sizes and genotypes from ranged reads of
`*_series_matrix.txt.gz` headers (`!Sample_characteristics_ch1`) and SOFT GSM records; deposit
formats and byte sizes from GEO FTP `suppl/` listings and `filelist.txt` (list and `HEAD` only);
cellxgene objects from the collection API plus h5ad header reads. No FASTQ, BAM, h5ad, `.rds`/
`Robj` or `RAW.tar` was downloaded.

**What was fetched.** Only the 1.27 MB `GSE270788_metadata.csv.gz` cell table. Every other size in
this note is a listing size, not a transfer.

**UNVERIFIED, and why it matters.**

- **No GEO series carries a licence field.** "GEO public" in the tables means NCBI GEO access
  terms; a per-paper CC-BY statement is plausible but was not verified. Confirm before
  redistributing any matrix inside the repo.
- **The Reichart ARVC/genotype column is known only from its `obs` schema**
  (`Primary.Genetic.Diagnosis`) — the ARVC (8) and NCCM (1) donor counts come from `obs`
  cross-tabs, not from a GEO series matrix.
- **Several CM-label sources live inside large objects** and were not opened:
  GSE298023's 4.0 GB `.rds`, GSE183852's 2.1–12 GB Seurat `Robj`, GSE185289's 2.5 GB CM Seurat
  object. For those rows the CM label is **UNVERIFIED** until the object is read.
- GSE101514 (per-arm n) and GSE136252 (CM-specific n) remain **UNVERIFIED**.
- GSE270788's AMI arm has **no per-sample time-since-infarct column**, so "acute" cannot be
  re-derived from the shipped metadata.

**Companion documents.** The dataset-level inventory is
[`stage2-stage3-gene-axis-datasets.md`](stage2-stage3-gene-axis-datasets.md). The fusion and
eligibility contract this survey feeds is ADR 0066
([`0066-stage2-corpus-view-and-gene-axis.md`](../adr/0066-stage2-corpus-view-and-gene-axis.md));
[`stage2-stage3-gene-axis-rank-aggregation.md`](stage2-stage3-gene-axis-rank-aggregation.md) is
the superseded ADR 0065-era aggregation evidence and is historical.
