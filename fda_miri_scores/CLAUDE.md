# CLAUDE.md — fda_miri_scores/

## Purpose

MIRI drug efficacy scoring from XPert-predicted cardiomyocyte transcriptome. Ranks 3,127 FDA-approved small-molecule drugs by their predicted ability to reverse the myocardial ischemia-reperfusion injury (MIRI) transcriptional disease signature. This is the final stage of the ARC drug discovery pipeline.

## Structure

```
fda_miri_scores/
├── xpert_fda-cm_deg.csv                         # XPert-predicted DEG: 3127 drugs × 978 landmark genes (logFC)
├── scoring.ipynb                                 # Scoring notebook: execution + visualization
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
│   └── MIRI_pathways.json                        # 16 MIRI-relevant KEGG pathways + direction + L1000 genes
├── scripts/                                      # Drug efficacy scoring pipeline
│   ├── __init__.py                               # Module docstring
│   ├── utils.py                                  # Shared: enrichment_score(), zscore_normalize(), data I/O
│   ├── score_wcs_rges.py                         # Score 1: WCS + RGES (gene-set based)
│   ├── score_zhang.py                            # Score 2: Zhang Score / ssCMap (rank-based)
│   ├── score_pathway_rges.py                     # Score 3: Pathway-level RGES (mechanistic)
│   ├── run_all.py                                # Main pipeline: all 3 methods → integrated ranking (4 ranks)
│   ├── generate_reference_genesets.py             # Generate all reference gene set JSONs from GEO data
│   └── build_union_and_rank.py                   # Build union gene set + ranked reference for Zhang
├── scores/                                       # Scoring output CSVs
│   ├── scores_wcs_rges.csv                       # WCS+RGES scores (all 3127 drugs)
│   ├── scores_zhang.csv                          # Zhang scores (all 3127 drugs)
│   ├── scores_pathway_rges.csv                   # Pathway-level scores (all 3127 drugs)
│   └── fda_miri_scores_integrated.csv            # Integrated 4-rank mean ranking
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
- **Reference**: Type C (MIRI_pathways.json, 16 pathways with direction)
- **Method**: Per-pathway RGES → efficacy profile. Disease-UP pathways: efficacy = −ES. Disease-DOWN pathways: efficacy = ES. Aggregated via mean_efficacy.
- **Key metric**: Higher mean_efficacy = better. n_effective = pathways with efficacy > 0.1.
- **Script**: `scripts/score_pathway_rges.py`

### Integration
- **Integrated rank**: Mean of 4 individual ranks (rank_WCS, rank_RGES, rank_Zhang, rank_pathway).
- **rank_pathway**: Descending (higher efficacy = better rank). All others: ascending.
- **Script**: `scripts/run_all.py`

## Key Files

### Input Data
- xpert_fda-cm_deg.csv — 3128 lines (header + 3127 drugs). Columns: drug_name, smiles, drug_idx, then 978 L1000 landmark genes. Values are predicted logFC from XPert ensemble. **Critical limitation**: 97.2% of values are negative (predominantly down-regulated). Z-score normalization (`--zscore`) strongly recommended.
- reference-genesets/MIRI_up-down.json — Union gene set (78 UP + 36 DOWN in L1000). Default gene set.
- reference-genesets/MIRI_rank.json — Ranked reference (107 genes from GSE160516) with signed ranks by |disease logFC|.
- reference-genesets/MIRI_pathways.json — 16 MIRI pathways (14 UP + 2 DOWN) with direction labels + L1000 genes.

### Output Scores
- scores_wcs_rges.csv — WCS, ES_up, ES_down, RGES, rank_WCS, rank_RGES, avg_rank
- scores_zhang.csv — Zhang_c, Zhang_C, Zhang_C_max, Zhang_p_approx, rank_Zhang
- scores_pathway_rges.csv — Per-pathway ES + efficacy, mean_efficacy, n_effective, frac_effective, rank_pathway
- fda_miri_scores_integrated.csv — Merged all-method scores + rank_integrated (mean of 4 ranks)

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

# Individual scoring
python scripts/score_wcs_rges.py --deg-csv xpert_fda-cm_deg.csv --gene-set reference-genesets/MIRI_up-down.json --zscore
python scripts/score_zhang.py --deg-csv xpert_fda-cm_deg.csv --rank-ref reference-genesets/MIRI_rank.json --zscore
python scripts/score_pathway_rges.py --deg-csv xpert_fda-cm_deg.csv --pathways reference-genesets/MIRI_pathways.json --zscore

# Regenerate reference gene sets
python scripts/generate_reference_genesets.py
python scripts/build_union_and_rank.py
```

## Key Findings (2026-06-07)

### Reference Gene Sets
- GSE160516 (mouse) is the primary reference: 75 UP + 32 DOWN L1000 genes
- GSE108940 (rat): only 10 L1000 genes (array annotation limits mapping)
- GSE137293 (pig): 3 L1000 genes (scar vs border, no sham control)
- Union gene set (MIRI_up-down.json): 78 UP + 36 DOWN, majority-vote conflict resolution
- 16 MIRI pathways curated (KEGG + literature), 14 UP + 2 DOWN

### Scoring Performance (with z-score normalization)
- **WCS**: 315/3127 drugs have WCS ≠ 0. Range [−0.41, +0.40]. 5 strong candidates (WCS < −0.2 AND RGES < −0.2).
- **RGES**: Wider discrimination than WCS. Lower = better reversal.
- **Zhang Score**: 696 drugs with c < −0.3. Best discrimination power among all methods.
- **Pathway RGES**: Best interpretability. Top drugs: Cabozantinib (mean_efficacy=0.513, 16/16 pathways), Ponatinib (0.454, 14/16), Ripretinib (0.436, 15/16) — all multi-target TKIs.
- **Cross-method agreement**: Generally low (ρ: 0.01–0.27), each method captures different aspects of drug response.
- **Consensus**: 0 drugs in top-50 across all 4 rank components → complementary nature of different scoring approaches.

### Recommendations
1. Always use `--zscore` flag for drug discrimination
2. MIRI_up-down.json is the default gene set (union, best coverage)
3. MIRI_rank.json preferred for Zhang Score (preserves disease magnitude)
4. Pathway scoring provides the most biological insight into mechanism
5. Integrated ranking (mean of 4 ranks) combines complementary perspectives
6. Upstream improvement: increase XPert model's inter-drug DEG variance

## Notes

- xpert_fda-cm_deg.csv is a copy from xpert/predict/fda-cm/xpert_fda-cm_deg.csv.
- GEO DEG results: mouse/rat via limma, pig via scipy Welch t-test.
- Pig data caveat: GSE137293 contains scar + border samples only (no sham control).
- The 978 genes correspond to the L1000 landmark gene set.
- Drug indices (drug_idx column) reference the shared 11,104-drug embedding index scheme.
- All 3 scoring methods share utils.py functions (enrichment_score, data I/O, z-score normalization).
