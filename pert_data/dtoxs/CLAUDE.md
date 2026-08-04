# CLAUDE.md — dtoxs/

## Purpose

DToxS (LINCS dataset LDS-1593) is a human cardiomyocyte drug perturbation transcriptomic dataset. It serves as the finetune target for the XPert model, adapting it from generic L1000 transcriptional space to cardiomyocyte-specific drug response prediction.

## Data Source

- Platform: 3'-DGE RNA-seq (Illumina HiSeq 2500, single-ended, 46 bp read length)
- Cell type: Primary human cardiomyocytes from PromoCell (4 donor lines: A, B, D, E; ages 54-62; 2 male, 2 female)
- Treatment: 48-hour drug exposure
- Samples: 914 FASTQ samples; 956 well-level records (206 control + 750 treated) across 5 sequencing experiments (Exps 1, 2, 4, 5; Exp 5 contains Cor, Spec-Fixed, Spec-Var-1/2/3 sub-experiments)
- GEO accession: GSE146096, GSE146097
- LINCS dataset ID: LDS-1593
- Reference: DToxS_article.md (Nature Communications, 2020)

## Directory Structure

```
dtoxs/
├── Data/
│   ├── Level3/
│   │   ├── Sequencing Experiment Number 1/    # Design + read-count TSVs
│   │   ├── Sequencing Experiment Number 2/
│   │   ├── Sequencing Experiment Number 4/
│   │   └── Sequencing Experiment Number 5/
│   │       ├── Cor/                            # Core experiment
│   │       ├── Spec-Fixed/
│   │       ├── Spec-Var-1/
│   │       ├── Spec-Var-2/
│   │       └── Spec-Var-3/
│   └── PromoCell-DGE-RNAseq-Experiment-Design.DCIC.xlsx
├── Metadata/
│   ├── Small_Molecule_Metadata.txt             # 48 small molecules, SMILES/InChI/batch info
│   ├── Antibody_Metadata.txt                   # Antibody perturbations (e.g., TRS/Trastuzumab)
│   ├── Protein_Metadata.txt                    # Protein perturbations (e.g., IGF, TNF)
│   ├── Primary_Cell_Metadata.txt               # 4 donor cardiomyocyte lines
│   ├── ExperimentalMetadata.txt                # Assay parameters
│   ├── DatasetMetadata.txt
│   ├── MetadataPropertyMapping.txt
│   ├── UCSC-hg38-Gene-HGNC-Entrez-Map.tsv      # Gene symbol/ID mappings
│   └── small_mol_category.md                   # Drug structural classification (6 categories)
├── ProcessingPipeline/
│   └── ProcessingPipelineSpecifications.txt
├── DToxS_article.md                            # Original publication summary
├── DToxS_summary.md                            # Well-level metadata coverage analysis
├── ReadMe.txt                                  # Original GEO package README
└── CLAUDE.md                                   # This file
```

## Drug Categories

DToxS contains 36 small molecules classified into 6 structural clusters (from small_mol_category.md):

- Category 1: Quinazoline-based KIs (5 drugs, 91 samples) — AFA, ERL, GEF, LAP, VAN
- Category 2: Pyrimidine-core KIs (10 drugs, 166 samples) — CER, DAB, DAS, IMA, NIL, PAZ, PON, RUX, TOF, TRA
- Category 3: Other heterocyclic KIs (8 drugs, 151 samples) — AXI, BOS, CAB, CRI, REG, SOR, SUN, VEM
- Category 4: Steroid-related (2 drugs, 30 samples) — EST, URS
- Category 5: Neuro/psychoactive drugs (5 drugs, 26 samples) — DOM, DTP, FLE, LOP, PAR
- Category 6: Other diverse therapeutics (6 drugs, 57 samples) — ALE, AMI, CYC, ENT, MIL, OLM

Total: 6 categories, 36 drugs, 521 small-molecule samples. Additionally, 750 treated wells include antibody (TRS) and protein (IGF, TNF) perturbations which lack SMILES.

This classification is used for the cold-cate (leave-one-category-out) cross-validation split in XPert finetuning.

## Processing Pipeline

The raw Level3 read counts are processed into miri_reposition/xpert-fda-v1/finetune/data/dtoxs.h5ad via these steps:

1. Map DTOXS read-count gene names to L1000 978 landmark genes using the UCSC-hg38 gene mapping table.
2. Compute log2(CPM + 1) per well.
3. Pair each treatment well with its corresponding control well, producing x_pert and x_ctl matrices.
4. Compute deg = x_pert - x_ctl (log fold change).

No executable code resides in this directory. The processing code lives in the xpert project.

## Relationship to XPert

- Input: Raw Level3 read counts are processed into miri_reposition/xpert-fda-v1/finetune/data/dtoxs.h5ad.
- Usage: Finetune target dataset for the XPert model (dual-branch Transformer).
- Data space: RMA-normalized log2 expression, non-negative, range approximately [0, 15].
- Compatibility with L1000: Both use log2 expression space. L1000 uses Q2NORM quantile-normalized log2; DTOXS uses RMA-normalized log2. No z-score conversion is applied to either dataset.
- Finetune strategies:
  - cold-drug: 10-fold leave-one-drug-out
  - cold-cate: 6-fold leave-one-category-out (using the 6 drug categories above)
  - dtoxs: single fold all-train (no validation/test, for final prediction model)
- Configs: dtoxs.yaml, cdsdb.yaml, panacea.yaml
- Split columns in h5ad: cold-drug_fold-1..10, cold-cate_fold-1..6, dtoxs_full

## Key Notes

- Gene mapping: DTOXS gene names differ from L1000 landmark genes in some cases due to gene symbol updates. A mapping exists for 27 genes (18 renamed + 9 direct matches). Examples: ELP1 was IKBKAP, COQ8A was ADCK3.
- Control states: CTRL, CTRL-DMSO, CTRL-EtOH (206 wells total across 5 experiments).
- Perturbation types: small molecules (672 wells with valid SMILES), antibodies (TRS, 19 wells), proteins (IGF 16 wells, TNF 14 wells), 29 wells unmatched.
- Cell lines: 4 primary human cardiomyocyte donors (A/B/D/E), all healthy Caucasians, ages 54-62, 2 male/2 female.
- 200 TSV files under Data/Level3/ (design + read-count files per sub-experiment).
- Data retrieved from DToxS center on 2020-10-17.
- dtoxs_full obs column is generated by miri_reposition/xpert-fda-v1/finetune/create_dtoxs_splits.py for final prediction model training.
