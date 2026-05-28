# CLAUDE.md — dtoxs/

## Overview

DToxS (LINCS dataset LDS-1593) — 人源心肌细胞药物扰动组数据集。用于 XPert 模型微调，使其适应心肌细胞环境下的药物转录组效应预测。

## Data Source

- **Platform**: 3'-DGE RNA-seq (Digital Gene Expression)
- **Cell type**: Human cardiomyocytes (PromoCell)
- **Treatment**: 48h drug exposure
- **Samples**: 914 samples (956 well-level records across 4 experiments)
- **GEO**: GSE146096, GSE146097
- **Reference**: [DToxS_article.md](DToxS_article.md)

## Directory Structure

```
dtoxs/
├── Data/Level3/                    # Raw read counts (TSV), organized by experiment
│   ├── Sequencing Experiment Number 1/
│   ├── Sequencing Experiment Number 2/
│   ├── Sequencing Experiment Number 4/
│   └── Sequencing Experiment Number 5/Cor/
├── Metadata/                       # Drug metadata, gene annotations
├── ProcessingPipeline/             # Original processing pipeline documentation
├── DToxS_article.md               # Original publication
├── DToxS_summary.md               # Dataset summary
└── ReadMe.txt                      # Original README from GEO
```

## Data Processing for XPert

The raw read counts are processed into `xpert/finetune/data/dtoxs.h5ad`:

1. Map DTOXS read-count gene names → L1000 978 landmark genes
2. Compute `log2(CPM + 1)` per well
3. Pair treatment wells with control wells → x_pert and x_ctl matrices
4. `deg = x_pert - x_ctl` (log fold change)

## Relationship to XPert

- **Input**: Raw Level3 read counts → processed into `xpert/finetune/data/dtoxs.h5ad`
- **Usage**: Finetune target dataset for XPert model
- **Data space**: RMA-normalized log2 expression (non-negative, range [0, ~15])
- **Compatible with L1000**: Both use log2 expression space (L1000: Q2NORM, DTOXS: RMA)

## Key Notes

- DTOXS gene names may differ from L1000 landmark genes due to gene symbol updates between projects. A mapping exists for 18 renamed genes (e.g., ELP1 formerly IKBKAP, COQ8A formerly ADCK3).
- Control states: CTRL, CTRL-DMSO, CTRL-EtOH
- No executable code in this directory — pure data and documentation
