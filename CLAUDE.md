# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a multi-project machine learning pipeline for drug-induced transcriptional perturbation prediction:

```
dtoxs/    →  Raw cardiotoxicity transcriptomic data (DToxS LINCS dataset LDS-1593)
kpgt/     →  Molecular graph transformer (KPGT) — extracts molecular features from SMILES
hg/       →  Heterogeneous GNN — generates biological drug embeddings via contrastive learning
xpert/    →  Main prediction model (XPert) — dual-branch transformer, consumes features from kpgt + hg
```

**Pipeline flow**: Drug SMILES → KPGT (molecular feats) + HG (biological embeddings) + UniMol (3D feats) → XPert (predicts 978-landmark-gene perturbation)

### Sub-project overviews

- **xpert/** — The active development project with detailed [CLAUDE.md](xpert/CLAUDE.md). Dual-branch transformer (~25M params) predicting post-drug transcriptional profiles. Has full pretrain/finetune/evaluate workflow.
- **kpgt/** — Molecular graph model (LiGhT architecture, DGL-based). Python 3.7, PyTorch 1.10, DGL-cu113. Conda env: `KPGT`. Workflow: pretrain on molecular graphs → extract features → downstream finetune.
- **hg/** — HeteroGNN (PyTorch Geometric) trained with InfoNCE contrastive loss on drug-gene-disease heterogeneous graph (PPI + DTI + DDS edges). Output: `(8981, 256)` drug embeddings consumed by XPert.
- **dtoxs/** — Pure data: 914 human cardiomyocyte 3'-DGE RNA-seq samples (956 well-level records), with metadata and read-count matrices. No executable code.

## Commands

### xpert/

```bash
# CLI training (from xpert/)
python -m src.train --mode train --config pretrain/config/l1000 --dataset l1000_mdmt_full --device cuda:0

# Shell scripts (env-var-driven)
bash pretrain/run.sh [dataset] [nfold] [config]
bash finetune/run.sh <pretrained_path> [dataset] [nfold] [config]

# Interactive pipeline
jupyter notebook run.ipynb          # 3-step: pretrain → finetune → evaluate
jupyter notebook metrics/evaluate.ipynb   # Evaluation only

# Config system: all paths in YAML resolved relative to xpert/ project root
# Outputs: pretrain/output/<dataset>/<timestamp>/  or  finetune/output/<dataset>-<fold>/<timestamp>/
```

### kpgt/

```bash
# Activate conda environment first: conda activate KPGT
cd kpgt/

# SMILES → KPGT pickle (feature extraction pipeline)
bash create_pickle.sh

# Pretrain molecular graph model
python scripts/train_kpgt.py --data_path <path> --save_path <path> --n_steps <N> --config base

# Extract features from pretrained model
python scripts/extract_features.py --config base --model_path <path> --data_path <path> --dataset <name>

# Finetune on downstream task
python scripts/finetune.py --config base --model_path <path> --data_path <path> --dataset <name> --dataset_type <classification|regression> --metric <name> --split <name>
```

### hg/

```bash
# Train HeteroGNN (from hg/)
bash run.sh                                    # default params
bash run.sh --device cuda:0 --epochs 500       # with overrides

# Or use the notebook: jupyter notebook run.ipynb

# Output: embeddings/HG_drug_embeddings.npy (8981, 256) — consumed by XPert
```

## Key Architectural Notes

- **XPert architecture**: See [xpert/CLAUDE.md](xpert/CLAUDE.md) for full details (model, config, data pipeline, key changes from history).
- **KPGT model**: LiGhT (Light Graph Transformer) — 12-layer attention over molecular graphs with virtual nodes. Uses DGL `HeteroGraphConv`-style message passing. Pretrained with path-based encoding and molecular descriptor reconstruction.
- **HG model**: HeteroGNN with 3 HeteroConv layers, neighbor sampling, InfoNCE contrastive loss over drug-gene edges. Uses `torch_geometric`.
- **Cross-project dependency**: KPGT features (`all_drugs_idx2KPGT.npy`) and HG embeddings (`all_drugs_idx2HG.npy`) are consumed by XPert as static drug representations in `drug_embeddings/`.
- **Environment isolation**: kpgt uses a dedicated conda env (`KPGT`, Python 3.7), while xpert and hg use system/modern Python (3.10+). They do not share a runtime — features are exchanged as pre-computed `.npy` files.
