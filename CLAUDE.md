# CLAUDE.md

## Project Goal

利用人工智能与生物信息学方法，预测可用于治疗心肌缺血再灌注损伤（MIRI）的药物。核心工作流：从小分子化学结构出发，预测其在心肌细胞上的转录组扰动效应（978个landmark基因）。

## Repository Structure

```
arc/
├── xpert/       → 药物转录组扰动预测模型（XPert），双分支Transformer，经L1000预训练+DTOXS微调
├── dtoxs/       → 心肌细胞药物扰动组数据集（DToxS LINCS LDS-1593），用于XPert模型微调
├── kpgt/        → 分子图Transformer（KPGT/LiGhT），将药物SMILES转化为分子图嵌入向量
├── hg/          → 异质图神经网络（HeteroGNN），生成药物的生物学嵌入向量（对比学习）
└── references/  → 虚拟细胞/药物预测相关参考文献（与当前工作流无关）
```

**Pipeline flow**: Drug SMILES → KPGT (molecular feats) + HG (biological embeddings) + UniMol (3D feats) → XPert (predicts 978-landmark-gene perturbation)

## Commands

### xpert/

```bash
cd xpert/
# CLI training
python -m src.train --mode train --config pretrain/config/l1000 --dataset l1000_mdmt_full --device cuda:0

# Shell scripts (env-var-driven)
bash pretrain/run.sh [dataset] [split_method] [config]
bash finetune/run.sh <pretrained_path> [dataset] [split_method] [config]

# FDA drug prediction pipeline
cd vcms/fda/
python query_fda_drugs.py           # Step 0: fetch FDA-approved drugs from ChEMBL
python step2_kpgt_prepare.py        # Step 1: prepare SMILES for KPGT
python step2_generate_embeddings.py # Step 2: generate UniMol + Morgan features
python gen_hg.py                    # Step 2b: generate HG biological embeddings
python step3_merge.py               # Step 3: merge new drug features into embedding arrays
python step4_ensemble_predict.py    # Step 4: 7-fold ensemble prediction
python postprocess_deg_rank.py      # Step 5: output deg + gene ranks per drug
```

### kpgt/

```bash
conda activate KPGT                # Python 3.7, DGL-cu113
cd kpgt/
bash create_pickle.sh               # SMILES → KPGT features
python scripts/train_kpgt.py --config base --data_path <path> --save_path <path>
python scripts/extract_features.py --config base --model_path <path> --data_path <path>
```

### hg/

```bash
cd hg/
bash run.sh                         # Train HeteroGNN
bash run.sh --device cuda:0 --epochs 500
# Output: model/<timestamp>/global_fold_early_stop.pth
```

## Cross-Project Dependencies

- kpgt → xpert: `xpert/drug_embeddings/all_drugs_idx2KPGT.npy` (molecular features)
- hg → xpert: `xpert/drug_embeddings/all_drugs_idx2HG.npy` (biological embeddings)
- kpgt uses conda env `KPGT` (Python 3.7); xpert and hg use conda env `xpert` (Python 3.13)
- Features are exchanged as pre-computed `.npy` files — no shared runtime

## Key Architectural Notes

- **XPert**: ~25M params, hidden=256, 8 attention heads. Dual-branch (pert + ctl) with cross-attention drug-cell fusion. See [xpert/CLAUDE.md](xpert/CLAUDE.md).
- **KPGT**: LiGhT architecture, 12-layer attention over molecular graphs with virtual nodes, DGL-based.
- **HG**: HeteroGNN, 3 HeteroConv layers, neighbor sampling, InfoNCE contrastive loss. PPI + DTI + DDS edges.
- **DTOXS data**: 914 human cardiomyocyte 3'-DGE RNA-seq samples, 48h drug treatment, RMA-normalized log2 expression.
- **L1000 data**: LINCS Phase I/II, Level 3 Q2NORM (quantile-normalized log2 expression), 978 landmark genes, non-negative.
- Both pretrain and finetune data are in log2 expression space — no z-score conversion needed.
