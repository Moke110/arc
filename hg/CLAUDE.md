# CLAUDE.md — hg/

## Overview

Heterogeneous Graph Neural Network (HeteroGNN) — 基于异质图的药物嵌入模型。在药物-基因-疾病异质图上通过对比学习训练，生成药物的生物学嵌入向量，作为 XPert 的药物特征输入之一。

## Environment

Uses conda env `xpert` (Python 3.13, PyTorch 2.12, PyTorch Geometric 2.7).

## Directory Structure

```
hg/
├── src/                              # Source code
│   ├── model_hg.py                   # HeteroGNN model definition
│   ├── pretrain_hg.py                # Training loop with InfoNCE loss
│   ├── utils_hg.py                   # Graph data loading, neighbor sampling
│   └── __init__.py                   # Package init
├── data/                             # Graph data (static .npy files)
│   ├── all_drug_node_feat_8981.npy           # Drug node features (8981 drugs)
│   ├── all_gene_node_feat_19392.npy          # Gene node features (19392 genes)
│   ├── PPI_edge_index_all_901260_pairs_weighted.npy  # Protein-protein interactions
│   ├── DTI_edge_index_all_12890_pairs.npy           # Drug-target interactions
│   ├── DDS_edge_index_all_287834_pairs_weighted.npy # Drug-drug similarities
│   └── 978_gene_node_idx.npy                  # L1000 landmark gene indices
├── model/                            # Trained model checkpoints
│   └── 20260518_102121/              # Latest training run
├── embeddings/                       # Output drug embeddings
├── run.sh                            # Training launcher (env-var-driven)
├── run.ipynb                         # Interactive training notebook
└── scratch/                          # Temporary/experimental files
```

## Key Workflows

### Training

```bash
cd hg/
bash run.sh                                    # Default params
bash run.sh --device cuda:0 --epochs 500 --lr 0.001
# Output: model/<timestamp>/global_fold_early_stop.pth
```

### Inference (generate embeddings for new drugs)

For FDA drugs, the inference pipeline is in `xpert/vcms/fda/gen_hg.py`:
1. Learns a ridge regression projection from KPGT (2304-dim) → HG input space (512-dim)
2. Runs HeteroGNN forward pass to get (256-dim) biological embeddings
3. Outputs: `fda_HG.npy`

## Model Architecture

- **HeteroGNN**: 3 `HeteroConv` layers with neighbor sampling
- **Node types**: Drug, Gene, Disease
- **Edge types**: PPI (gene-gene), DTI (drug-gene), DDS (drug-drug)
- **Loss**: InfoNCE contrastive loss over drug-gene edges
- **Hidden dim**: 256
- **Output**: (N_drugs, 256) drug embeddings

## Graph Data

The heterogeneous graph is pre-built from public databases:
- **PPI**: 901,260 weighted protein-protein interaction pairs
- **DTI**: 12,890 drug-target interaction pairs
- **DDS**: 287,834 weighted drug-drug similarity pairs
- **Drug nodes**: 8,981 drugs with initial features
- **Gene nodes**: 19,392 genes with initial features; 978 are L1000 landmarks

## Relationship to XPert

HG embeddings are consumed by XPert as `xpert/drug_embeddings/all_drugs_idx2HG.npy` (8981×256). For inference on new drugs, `gen_hg.py` generates embeddings via a learned linear projection from KPGT space → HG input space.
