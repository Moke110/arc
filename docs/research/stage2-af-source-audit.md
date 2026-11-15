# Stage 2 gene axis — the AF source audit (why AF stays out)

**Status:** research note (read-only audit, 2026-11-15). Not an ADR; authorises no download into
an external item. Reference: ADR request of 2026-11-14 to re-admit `AF` with `GSE255612` as its
CM DEG source after the "AF data search" session, ADR 0066 *Source eligibility*, ADR 0068
*Degenerate sources*, ADR 0069 *Remove the AF category*, and the audit decision ADR 0070.

**Superseded by ADR 0071 (2026-11-15), for its decision only.** The audit's *conclusion* — no
clinical AF dataset carries AF-specific ordering strong enough to make a fusion source — is
what ADR 0071 overturns: `AF` is published with three CM-specific tier-1 sources and a tier-3
atrial-bulk descent, and its honest `DEG2000` is the single gene `SLC9A7`. The audit's
*measurements* are still the source of record for the items it examined (they are cited by
ADR 0071 decision 7 and by the `AF` category catalog's exclusion reasons), so it is kept
rather than deleted. Read it as the input that ADR 0071 overrode, never as a live argument for
keeping `AF` out.

## Question

ADR 0069 removed the `AF` category because its only dataset, `GSE255612` (`gse255612__af`), was
non-informative: every axis-measured effect score `|log2FC|·−log10(padj)` was exactly 0. The user
then asked to re-admit AF with `GSE255612` confirmed as a clinical-sample AF source. Before
re-admitting, two questions had to be settled on data rather than provenance:

1. Were all AF sources found by the search actually excluded, and on what evidence?
2. Under the user's rule — **the data must be able to resolve cardiomyocytes specifically** — does
   any AF dataset carry an AF-specific ordering that could justify a fusion source?

## Method

- **Candidate list.** The AF-data-search session verified 117 GEO series individually
  (`tmp`/scratch artifacts, untracked) plus a hiPSC/hESC-CM model inventory. Every accession was
  cross-checked against the tracked repo (external items, ADRs, research notes, the AF catalog's
  `rejected[]`).
- **CM-specificity filter (user rule).** Whole-tissue bulk (including GSE269882, GSE128188,
  GSE226282-4, GSE161450, GSE127856, GSE79768 and the controlled-access EGA bulk cohort),
  non-CM cell types (fibroblasts, CD45+/immune, endothelial, non-CM fractions), blood/PBMC,
  microarray, methylation/spatial/ChIP and unlabelled atlases are out of scope **whatever their
  DEG count**: they cannot resolve a cardiomyocyte-specific transcriptome.
- **Actual DEG counts.** For each remaining candidate the processed payload was downloaded to
  scratch and re-analysed at donor level (limma-voom; edgeR QL and DESeq2 as sensitivity), BH
  padj. No item tree was modified.
- **Label-permutation null.** For every candidate DE table, the top-`N` (`N` = 2,000/5,000)
  overlap with the released axes was compared with the null distribution obtained by permuting
  the sample group labels (100 permutations, same model refit). A uniform-random gene baseline is
  *not* used as the test: it flatters any dataset whose measurable genes happen to be large-effect
  genes.
- **Publication cross-check.** `GSE255612`'s own CM DEG table (Supplementary Data 4 of PMID
  39562555) was downloaded and read directly.

## Results — clinical, CM-specific

| source | design | n | CM access | actual CM DEGs (padj<0.05) | label-permutation null | verdict |
|---|---|---|---|---|---|---|
| **GSE255612** | human LA 10x 3′ v3 snRNA-seq; 179,697 nuclei, 34 donors — *GEO `overall_design` says 19/19; the metadata and the paper give 18 AF / 16 control (19+17 recruited, 2 failed QC)* | 18 AF / 16 control | author `cell_type == "Cardiomyocytes"`, 32,438 nuclei | **0** in this repo's pipeline (min padj 0.7359; all 18,092 axis-measured effect scores 0); **2** in the paper (LINC01257 +1.44, ATRNL1 +1.58, adj p 0.0174); count-model re-derivation: limma-voom **0**, edgeR QL **2**, DESeq2 **6**, union not credible AF biology (PKP1, MT-CO1, ACVRL1, ATRNL1, P3H2, LINC01257) and NPPA/MYH7/TTN/MYH6/NPPB all null | within-axis top-2000 **p=0.634** (320 vs 329.4±34.1), top-5000 **p=0.347**; genome-wide p=0.515 / 0.277 | **no AF-specific signal** |
| **GSE138252** | bulk RNA-seq of FANS-sorted PCM1+ CM nuclei (HiSeq 2500); LA/RA/LV/RV | LA 5 AF+HF / 5 NF; RA 5/3 (RA sex-confounded, unused) | CM-specific nuclear isolation, verified four ways (series summary, every GSM extract protocol, `molecule subtype`, paper) | LA 5v5: limma-voom **2**, edgeR QL **6** (LRRC10 −2.51, SLC9A7 +2.26, HBB +5.57, CHST7 +2.19, HOTAIRM1 −3.16, TRAC −3.56) | within-axis **p=0.574** (341 vs 336.0±35.2) and **p=0.901** (1,521 vs 1,630.3±112.9, below null); hits track HF not AF (AF-vs-CTL vs HF-only log2FC r=0.577; sex-adjusted FDR<0.05 = 0) | **no AF-specific signal**; AF+end-stage-HF confound is not de-confoundable (HF-only arm is n=3 and batch-separated) |
| **GSE238242** | human LAA 10x multiome (snRNA arm), persistent AF vs SR | 4 SR / 3 AF donors, 3,049 CM nuclei | author `cell_type == "CM"` | **0** (min padj 0.217; voom same; filter/trend sensitivity all 0) | not run (null DE) | no signal at n=4/3 |
| **GSE271229** | human LAA snRNA-seq, 118,275 nuclei | 3 no-AF / 3 paroxysmal / 3 persistent | **none** — the deposit ships only per-sample CellRanger barcodes/features/matrix; no cell-type annotation of any kind | not testable without inventing a gate | — | **ineligible** |
| dbGaP phs003644 | left-atrial CM RNA-seq, post-op AF | — | CM | — | — | controlled access |

`GSE255612`'s GEO series is internally inconsistent (`overall_design` 19/19 vs 18/16); the
per-nucleus metadata and the paper's own tables both give 18/16, which is what the audit uses.

## Results — hiPSC/hESC-CM model datasets (tier 2)

These are CM by construction, so they pass the user's CM rule, but they are in-vitro perturbation
models, not clinical AF. Counts are log2(CPM+1)+limma (voom in parentheses) at padj<0.05.

| dataset | contrast | n (case/control) | DEG count | note |
|---|---|---|---|---|
| GSE175944 | PITX2−/− vs control atrial EHT | 3/3 | **586** (501; batch-adjusted 440) | credible atrial model (ECM programme); PITX2 log2FC −6.78 = knockout confirmed |
| GSE255991 | ATRNL1 siRNA vs scramble / ATRNL1 overexpression | 3/6 (pooled); 3/3 | **1,963** (2,455); OE **7,783** | credible perturbation model, but it tests a gene knockdown, not AF disease |
| GSE307021 | TBX5-degron +dTAG 24 h vs DMSO | 3/3 | **302** (978) | credible atrial model (loss of atrial identity) |
| GSE221262 | FAM13B knockdown vs control | 3/3 | 458 | downgraded: normalised matrix only (no counts); iCell CM are ventricular-like |
| GSE271839 | LINC01629 CRISPRi vs NTC | 3/3 | **0** | knockdown real, no informative DEGs |
| GSE154084 | SCN5A E428K vs control/corrected | 2/2 | not usable | differentiation-state/line confound (pluripotency up, sarcomere down); matrix columns unlabelled |
| GSE312917 | TTN T32756I vs WT | 3/3 | 0 (submitter table) | no counts deposited; WT and mutant are different cell lines, not isogenic |
| GSE128908, GSE293813, GSE285169 | MYL4, PITX2c, TBX5 | 1–2 | not testable | n=1 (MYL4 has one biological replicate; PITX2c one per arm); TBX5 n=2 and a 1 GB deposit, not quantified |

`GSE234141` (TTNtv hiPSC-CM) is already declared inside `genetic_cardiomyopathy`; it is not an AF
source.

## Conclusions

1. **No CM-specific clinical AF dataset carries an AF-specific ordering.** The two largest
   (`GSE255612`, 18/16; `GSE138252`, 5/5) are indistinguishable from their label-permutation
   nulls; the two smaller ones are null at n=4/3 or not resolvable at all; the largest clinical
   deposit (dbGaP phs003644) is controlled-access.
2. **`GSE255612`'s flatness is a property of the data, not of the normalisation.** The paper's own
   CM DEG list has two genes and the authors state power was insufficient; a count-based
   re-derivation adds a handful of non-credible hits but no axis-relevant signal.
3. **Only tier-2 in-vitro models produce large DEG lists**, and they are perturbation models of AF
   GWAS genes (PITX2, ATRNL1, TBX5), not clinical AF. Declaring them would misdescribe the
   category.
4. **ADR 0068's informativeness rule is necessary but not sufficient.** A count-based `GSE255612`
   table would pass it (14,460 distinct non-zero effect scores) while carrying no AF-specific
   signal. ADR 0070 adds the label-permutation null as the admission standard for any new
   category dataset.

## Consequences for the products

- `AF` stays out of `01_bases` (ADR 0069 stands; ADR 0070 records the audit) and is a declared gap
  like `HHD`/`PPCM`.
- The gene-axis rule document gains a **source dataset information table** at its start (see
  ADR 0070); its second block lists the audited-and-excluded AF sources with the per-dataset DEG
  counts above, so the exclusion is visible in the document, not only in this note.
- DEG counts in that table use the **published `catalog.json` `n_significant`** (one row per axis
  gene, best contrast, padj<0.05) — the same number the release publishes. Item-level row counts
  are a different quantity (e.g. `GSE234141`'s four contrasts give 331 significant rows against
  the released 121 genes) and must not be mixed into the table.

## Audit artifacts (scratch, untracked; regenerate with the listed commands)

- `.scratch/af-public-data/` — the search session's inventory (`REPORT.md`, `clinical.md`,
  `models.md`, `verified_geo_review.tsv`).
- `.scratch/af-readmit-probe/REPORT.md` — `GSE255612` count-based re-derivation + permutation null.
- `.scratch/af-source-audit/GSE138252/{summary.md,permutation.md}` — `GSE138252` DE + null.
- `.scratch/af-source-audit/GSE238242/REPORT.md`, `.scratch/af-source-audit/GSE271229/REPORT.md`.
- `.scratch/af-source-audit/models/AF_MODEL_DEG_REPORT.md` — model DEG counts.
- `.scratch/verify-paper/deg.xlsx` — the paper's Supplementary Data 4 (2 CM DEGs).

## See

- [0066](../adr/0066-stage2-corpus-view-and-gene-axis.md) — source eligibility
- [0068](../adr/0068-fusion-amendment-degenerate-sources-and-rank-norm.md) — the non-informative rule
- [0069](../adr/0069-remove-af-category-from-01-bases.md) — the removal this audit confirms
- [0070](../adr/0070-af-stays-out-of-the-gene-axis.md) — the decision record
- [stage2-gene-axis-cm-source-survey.md](stage2-gene-axis-cm-source-survey.md) — the earlier CM-source survey
