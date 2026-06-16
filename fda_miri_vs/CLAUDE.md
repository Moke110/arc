# CLAUDE.md — fda_miri_scores/

## Purpose

MIRI drug efficacy scoring from XPert-predicted cardiomyocyte transcriptome. Ranks 3,127 FDA-approved small-molecule drugs by their predicted ability to reverse the myocardial ischemia-reperfusion injury (MIRI) transcriptional disease signature. This is the final stage of the ARC drug discovery pipeline.

## Structure

```
fda_miri_scores/
├── xpert_fda-cm_deg.csv                         # XPert-predicted DEG: 3127 drugs × 978 landmark genes (logFC)
├── fda_screen.ipynb                              # Drug screening notebook: scoring → filtering → ranking
├── scores/fda_screened.csv                       # Final ranked MIRI drug candidates
├── CLAUDE.md                                     # This file
├── README.md                                     # Human-readable Chinese documentation
├── geo/                                          # MIRI experimental transcriptome data + preprocessing
│   ├── GSE108940_WT_deg_results.csv              # Rat MIRI model DEG (28,037 genes, limma)
│   ├── GSE160516_deg_results.csv                 # Mouse MIRI model DEG (20,710 genes, limma)
│   ├── GSE160516_up_l1000.csv                    # MIRI UP genes → human L1000 (75 genes)
│   ├── GSE160516_down_l1000.csv                  # MIRI DOWN genes → human L1000 (32 genes)
│   ├── GSE108940_human_l1000.csv                 # Rat MIRI DEGs → human L1000 via homologene (10 genes)
│   ├── GSE137293_COUNTS.txt.gz                   # Pig I/R RNA-seq raw counts (scar + border, no sham)
│   ├── GSE137293_deg_results.csv                 # Pig scar vs border DEG (scipy t-test, 13,117 genes)
│   ├── GSE137293_human_l1000.csv                 # Pig DEGs → human L1000 (0 genes via homologene)
│   ├── run_deg_analysis.R                        # Mouse+rat DEG analysis (limma)
│   ├── map_to_human.R                            # Mouse→human ortholog mapping (homologene)
│   └── map_orthologs.R                           # Rat+pig→human ortholog mapping (homologene)
├── reference-genesets/                           # Generated JSON reference gene sets
│   ├── GSE160516_MIRI.json                       # Primary: mouse MIRI (75 UP + 32 DOWN, L1000)
│   ├── GSE108940_MIRI.json                       # Rat MIRI (4 UP + 6 DOWN, L1000)
│   ├── GSE137293_scar_vs_border.json             # Pig scar vs border (2 UP + 1 DOWN, L1000)
│   ├── MIRI_up-down.json                         # Union gene set (78 UP + 36 DOWN, all 3 species)
│   ├── MIRI_rank.json                            # Ranked reference for Zhang Score (107 genes, GSE160516)
│   └── MIRI_pathways.json                        # 40 MIRI pathways (20 UP + 20 DOWN), programmatically built
├── scripts/                                      # Drug efficacy scoring pipeline
│   ├── __init__.py                               # Module docstring
│   ├── utils.py                                  # Shared: enrichment_score(), zscore_normalize(), data I/O, coverage filter
│   ├── score_wcs_rges.py                         # Score 1: WCS + RGES (gene-set based)
│   ├── score_zhang.py                            # Score 2: Zhang Score / ssCMap (rank-based)
│   ├── score_pathway_rges.py                     # Score 3: Pathway-level RGES (mechanistic)
│   ├── run_all.py                                # Main pipeline: all 3 methods → integrated ranking (4 ranks)
│   ├── generate_reference_genesets.py             # Generate per-dataset gene set JSONs from GEO data
│   ├── build_union_and_rank.py                   # Build union gene set + ranked reference for Zhang
│   └── rebuild_pathways.py                       # Rebuild MIRI_pathways.json from GEO GSEA + literature
├── .claude/                                      # Intermediate files (cache, enrichment results, logs)
├── scores/                                       # Scoring output CSVs
│   ├── scores_wcs_rges.csv                       # WCS+RGES scores + WCS_zs, RGES_zs (3127 drugs)
│   ├── scores_zhang.csv                          # Zhang scores + Zhang_zs (3127 drugs)
│   ├── scores_pathway_rges.csv                   # Pathway-level scores + pathway_reversal_zs (3127 drugs)
│   ├── fda_miri_scores_integrated.csv            # Integrated 4-rank mean ranking
│   └── fda_filtered_scores.csv                   # Drugs passing WCS<0 & all critical pathways>0
└── plots/                                        # Scoring visualization PNGs
    ├── score_distributions.png                   # 5 score distributions + agreement histogram
    ├── top20_integrated.png                      # Top 20 drug ranking + 5-method heatmap
    ├── method_agreement.png                      # 6×6 Spearman ρ heatmap
    ├── method_pairwise_scatter.png               # 8 pairwise score scatter plots
    ├── pathway_heatmap_top10.png                 # Top 10 drug pathway reversal heatmap
    ├── top50_consensus.png                       # Top 50 consensus (Venn + bar chart)
    └── top_vs_bottom.png                         # Top 50 vs Bottom 50 score comparison
```

## Three Scoring Pipelines

### Pipeline 1: WCS + RGES (gene-set based)
- **Reference**: Type A (up/down gene sets, e.g. MIRI_up-down.json)
- **WCS**: Weighted Connectivity Score. (ES_up − ES_down)/2 when signs differ, else 0. Subramanian et al. 2017 *Cell*.
- **RGES**: Reverse Gene Expression Score. |ES_up| − |ES_down|. Always defined; Chen et al. 2017 *Nat. Commun.*.
- **Key metric**: Lower = better reversal. WCS < −0.2 AND RGES < −0.2 = strong candidate.
- **Script**: `scripts/score_wcs_rges.py`

### Pipeline 2: Zhang Score / ssCMap (rank-based)
- **Reference**: Type A (up/down gene sets) or Type B (MIRI_rank.json, preferred)
- **Method**: Signed rank product. Drug genes ranked by |logFC| → signed ranks R(g). Disease genes with pre-computed signed ranks s(g). c = Σ R(g)·s(g) / C_max. Zhang & Gant 2008 *BMC Bioinformatics*.
- **Key metric**: c ∈ [−1, 1], negative = reversal. c < −0.3 = strong candidate.
- **Script**: `scripts/score_zhang.py`

### Pipeline 3: Pathway-Level RGES (mechanistic)
- **Reference**: Type C (MIRI_pathways.json, 40 pathways: 20 UP + 20 DOWN)
- **Method**: Per-pathway RGES → efficacy profile. Disease-UP pathways: efficacy = −ES. Disease-DOWN pathways: efficacy = ES. Aggregated via mean_efficacy.
- **Key metric**: Higher mean_efficacy = better. n_effective = pathways with efficacy > 0.1.
- **Script**: `scripts/score_pathway_rges.py`
- **Rebuild script**: `scripts/rebuild_pathways.py` — programmatic pathway generation from GEO GSEA + literature

### Integration
- **Legacy integration**: Mean of 4 individual ranks (rank_WCS, rank_RGES, rank_Zhang, rank_pathway). Script: `scripts/run_all.py`.
- **Current screening pipeline** (fda_screen.ipynb):
  1. **Scoring**: 3 methods → 4 Z-score columns + 6 critical pathway efficacy Z-scores (e.g. `KEGG-Apoptosis_efficacy_zs`).
  2. **Filtering**: WCS < 0 AND ≥4/6 critical pathway efficacies > 0 → 199 drugs pass → `scores/fda_filtered_scores.csv`.
  3. **Ranking**: integrated_score = −0.3×RGES_zs − 0.1×Zhang_zs + 0.3×pathway_reversal_zs + 0.05×Σcritical_efficacy_zs → `scores/fda_screened.csv`.

## Key Files

### Input Data
- xpert_fda-cm_deg.csv — 3128 lines (header + 3127 drugs). Columns: drug_name, smiles, drug_idx, then 978 L1000 landmark genes. Values are predicted logFC from XPert ensemble. **Critical limitation**: 97.2% of values are negative (predominantly down-regulated). Z-score normalization (`--zscore`) strongly recommended.
- reference-genesets/MIRI_up-down.json — Union gene set (78 UP + 36 DOWN in L1000). Default gene set.
- reference-genesets/MIRI_rank.json — Ranked reference (107 genes from GSE160516) with signed ranks by |disease logFC|.
- reference-genesets/MIRI_pathways.json — 40 MIRI pathways (20 UP + 20 DOWN) with KEGG-/GO-/KEGG&GO- naming. Generated by `scripts/rebuild_pathways.py`.

### Output Scores
- scores_wcs_rges.csv — WCS, ES_up, ES_down, RGES, WCS_zs, RGES_zs, rank_WCS, rank_RGES, avg_rank
- scores_zhang.csv — Zhang_c, Zhang_C, Zhang_C_max, Zhang_p_approx, Zhang_zs, rank_Zhang
- scores_pathway_rges.csv — Per-pathway ES + efficacy, mean_efficacy, n_effective, pathway_reversal, pathway_reversal_zs, 6 critical pathway efficacy Z-scores, rank_pathway
- fda_miri_scores_integrated.csv — Merged all-method scores + rank_integrated (mean of 4 ranks)
- fda_filtered_scores.csv — Drugs passing WCS<0 AND ≥4/6 critical pathway efficacies > 0 (199 drugs); includes WCS_zs, RGES_zs, Zhang_zs, pathway_reversal_zs + per-pathway efficacy
- scores/fda_screened.csv — Final ranked drug candidates with integrated_score = −0.3×RGES_zs − 0.1×Zhang_zs + 0.3×pathway_reversal_zs + 0.05×Σcritical_efficacy_zs

## Dependencies

- xpert/predict/fda-cm/ — Source of xpert_fda-cm_deg.csv
- GEO public data — GSE108940, GSE160516, GSE137293
- Python: numpy, pandas, scipy, tqdm, matplotlib, seaborn
- R: homologene (for species ortholog mapping)

## Scoring Usage

```bash
# Full pipeline (all 3 methods, z-scored, outputs to scores/)
cd fda_miri_scores/
python scripts/run_all.py --zscore

# Drug screening notebook (score → filter → rank)
jupyter notebook fda_screen.ipynb

# Individual scoring
python scripts/score_wcs_rges.py --deg-csv xpert_fda-cm_deg.csv --gene-set reference-genesets/MIRI_up-down.json --zscore
python scripts/score_zhang.py --deg-csv xpert_fda-cm_deg.csv --rank-ref reference-genesets/MIRI_rank.json --zscore
python scripts/score_pathway_rges.py --deg-csv xpert_fda-cm_deg.csv --pathways reference-genesets/MIRI_pathways.json --zscore

# Regenerate reference gene sets
python scripts/generate_reference_genesets.py
python scripts/build_union_and_rank.py

# Rebuild pathway gene sets from GEO enrichment + literature
python scripts/rebuild_pathways.py [--target-up 30] [--target-down 30] [--skip-literature]
```

## Key Findings (2026-06-12)

### Reference Gene Sets
- GSE160516 (mouse) is the primary reference: 75 UP + 32 DOWN L1000 genes
- GSE108940 (rat): only 10 L1000 genes (array annotation limits mapping)
- GSE137293 (pig): 3 L1000 genes (scar vs border, no sham control)
- Union gene set (MIRI_up-down.json): 78 UP + 36 DOWN, majority-vote conflict resolution
- **MIRI_pathways.json rebuilt (2026-06-11)**: 40 pathways (20 UP + 20 DOWN), programmatically generated from GEO GSEA (GSE160516/GSE108940/GSE137293) against KEGG_2021_Human (320 pathways) and GO_Biological_Process_2023 (5,407 terms) via Enrichr. 8,619 candidate pathways evaluated, 1,660 passed dual filters (n_genes > 5 AND coverage ≥ 978/20000), 53 merged. Two-phase selection: literature-supported pathways prioritized, then GEO-filled to 20/direction. 20 literature-validated pathways (13 UP + 7 DOWN).

### Pathway Building Pipeline
- **Script**: `scripts/rebuild_pathways.py`
- **Sources**: GEO DEG enrichment (preranked GSEA using historical enrichment_score) + gpt-researcher literature survey
- **Ortholog mapping**: mygene + homologene (mouse/rat/pig → human Entrez → human symbol)
- **Filters**: (1) n_l1000 > 5; (2) n_l1000/n_total ≥ 978/20000 = 0.0489
- **Merging**: KEGG+GO pathways with identical L1000 gene content + same direction → KEGG&GO- prefix
- **Naming**: All pathways prefixed with KEGG- or GO- (or KEGG&GO- for merged)
- **Selection**: Two-phase prioritization. Phase 1: literature-supported pathways selected first (sorted by composite score). Phase 2: remaining slots filled by top-scored GEO-only pathways. Composite score = GEO ES magnitude × dataset weight + literature support + gene count bonus.
- **Caching**: Enrichr gene sets cached 1 week, ortholog mappings permanent in `.claude/`
- **utils.py additions**: `load_l1000_genes()`, `passes_coverage_filter()`; `load_pathway_gene_sets()` min_genes updated to 6

### Biological Themes in New Pathways
- **UP in MIRI** (20): Inflammatory/immune signaling (TNF, NF-κB, MAPK, chemokine, IL-1β production), apoptosis/necroptosis, hypoxia/oxidative stress response, wound healing, leukocyte adhesion/migration, T cell activation
- **DOWN in MIRI** (20): Cardiac muscle contraction, calcium signaling, adrenergic signaling, focal adhesion/ECM organization, actin cytoskeleton regulation, dilated cardiomyopathy, fatty acid degradation, peroxisome, glycolysis, pyruvate metabolism, mitochondrial organization, insulin signaling, TORC1 signaling

### Scoring Performance (with z-score normalization)
- **WCS**: 315/3127 drugs have WCS ≠ 0. Range [−0.41, +0.40].
- **RGES**: Wider discrimination than WCS. Lower = better reversal.
- **Zhang Score**: 696 drugs with c < −0.3. Best discrimination power among all methods.
- **Pathway RGES**: Best interpretability. pathway_reversal = fraction of 40 pathways with efficacy > 0.1.
- **Z-scores**: WCS_zs, RGES_zs, Zhang_zs, pathway_reversal_zs — standardized across drugs for cross-method comparison and integration.
- **Cross-method agreement**: Generally low (ρ: 0.01–0.27), each method captures different aspects of drug response.

### Screening Pipeline (fda_screen.ipynb)
- **Filtering**: WCS < 0 AND ≥4/6 critical pathway efficacies > 0 → 199/3127 drugs pass.
- **Ranking**: integrated_score = −0.3×RGES_zs − 0.1×Zhang_zs + 0.3×pathway_reversal_zs + 0.05×Σcritical_efficacy_zs. Sign convention: RGES_zs, Zhang_zs more negative = better → negated. pathway_reversal_zs more positive = better → positive weight.
- **Top candidates**: 1. Metocurine Iodide (integrated=+0.008), 2. Penicillin G (−0.047), 3. Nalidixic Acid (−0.096), 4. Leucovorin Calcium (−0.106).
- **Output**: `scores/fda_screened.csv` — 199 final ranked MIRI drug candidates.

### Recommendations
1. Always use `--zscore` flag for drug discrimination
2. MIRI_up-down.json is the default gene set (union, best coverage)
3. MIRI_rank.json preferred for Zhang Score (preserves disease magnitude)
4. Pathway scoring provides the most biological insight into mechanism
5. Screening pipeline (fda_screen.ipynb) combines 4 Z-score components with weighted integration
6. Upstream improvement: increase XPert model's inter-drug DEG variance

## Notes

- xpert_fda-cm_deg.csv is a copy from xpert/predict/fda-cm/xpert_fda-cm_deg.csv.
- GEO DEG results: mouse/rat via limma, pig via scipy Welch t-test.
- Pig data caveat: GSE137293 contains scar + border samples only (no sham control).
- The 978 genes correspond to the L1000 landmark gene set.
- Drug indices (drug_idx column) reference the shared 11,104-drug embedding index scheme.
- All 3 scoring methods share utils.py functions (enrichment_score, data I/O, z-score normalization).
- **gpt-researcher config**: Uses local embedding model at `.agents/skills/gpt-researcher/models/all-MiniLM-L6-v2/` (87.4 MB). LLM: DeepSeek (FAST_LLM/SMART_LLM/STRATEGIC_LLM in .env). HuggingFace mirror: HF_ENDPOINT=https://hf-mirror.com.
