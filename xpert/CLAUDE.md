# CLAUDE.md

XPert is a dual-branch transformer for predicting drug-induced transcriptional perturbation responses (978 landmark genes). It takes drug molecular features (UniMol + HG biological embeddings + dose/time tokens) and cell state (binned expression + PPI gene vectors), and outputs predicted perturbed expression, control expression, and expression change (x_deg).

> **历史项目文档**: 见 [.claude/xpert-history_summary.md](.claude/xpert-history_summary.md) — 原始代码仓库 file-by-file 详细总结

## Project Structure

```
xpert/
├── run.ipynb                  # 3-step pipeline launcher (pretrain → finetune → evaluate)
├── src/                       # All source code
│   ├── train.py               # CLI entry point, arg parsing, train/validate/infer loops
│   ├── model_xpert.py         # XPertNet + AttnEncoder (dual-branch transformer)
│   ├── model_utils.py         # Attention, Embeddings, Encoder building blocks (uses F.scaled_dot_product_attention)
│   ├── dataset.py             # MyDataset: AnnData → PyTorch tensors (lazy __getitem__ with pre-extracted numpy arrays)
│   ├── utils.py               # Data loading, loss functions, EarlyStopping, UniMol feat gen
│   ├── config.py              # YAML config loader, path resolver, output dir helper
│   ├── metrics.py             # Core metrics (MSE, PCC, Precision@K) + get_metrics_infer
│   ├── eval_metrics.py        # Extended metrics (R², MMD, Wasserstein) + get_metrics_new
│   ├── vis_gene_attn.py       # Gene attention heatmaps/boxplots/scatter
│   ├── vis_molecule.py        # Molecule attention visualization with RDKit
│   └── __init__.py            # Lazy imports
├── pretrain/                  # Step 1: pretrain
│   ├── config/l1000.yaml      # Pretrain config (batch 64, lr 0.004)
│   ├── data/                  # Input data (h5ad)
│   ├── output/                # Pretrained model weights
│   └── run.sh                 # Launch script (env-var-driven)
├── finetune/                  # Step 2: finetune
│   ├── config/cdsdb.yaml      # Clinical CDSDB config (batch 32)
│   ├── config/panacea.yaml    # PANACEA config (batch 32)
│   ├── data/                  # Target dataset h5ad files
│   ├── output/                # Finetuned model weights
│   └── run.sh                 # Launch script
├── metrics/                   # Step 3: evaluate
│   ├── evaluate.ipynb         # Interactive evaluation notebook
│   └── output/                # Evaluation results (CSV, npy, pkl)
├── drug_embeddings/           # Shared drug feature files
│   ├── all_drugs_unimol_arr.npy       # UniMol 3D molecular features (8981×122×514)
│   ├── all_drugs_idx2HG.npy           # Pre-trained HG drug embeddings (8981×256)
│   ├── PPI_gene_vector_128d.npy       # PPI gene network positional embeddings (978×128)
│   ├── all_drugs_idx2smi_8981.npy     # Drug index → SMILES mapping
│   └── all_drugs_idx2KPGT.npy / all_drugs_idx2morgan.npy  # Alternative drug features
└── metadata/                  # Gene/drug reference info
```

## Key Files to Read for Context

| When asked about... | Read these files |
|---|---|
| Model architecture | `src/model_xpert.py`, `src/model_utils.py` |
| Training loop / CLI args | `src/train.py` (arg_parse + main) |
| Data pipeline / h5ad format | `src/dataset.py`, `src/utils.py` (_load_data_and_drug_feat) |
| Config / hyperparameters | `pretrain/config/l1000.yaml`, `finetune/config/cdsdb.yaml`, `src/config.py` |
| Loss functions | `src/train.py` (train function), `src/utils.py` (mse_loss_ls_sum, pcc_loss_sum) |
| Evaluation metrics | `src/eval_metrics.py` (get_metrics_new), `src/metrics.py` (get_metrics_infer) |
| Visualization | `src/vis_gene_attn.py`, `src/vis_molecule.py` |
| Pipeline workflow | `run.ipynb`, `pretrain/run.sh`, `finetune/run.sh` |
| FDA inference pipeline | `vcms/fda/` — query FDA drugs → embeddings → ensemble predict → postprocess deg/rank |

## Entry Points

- **CLI**: `python -m src.train --mode train --config pretrain/config/l1000 --dataset l1000_mdmt_full ...`
- **Notebook**: `run.ipynb` (3-step launcher), `metrics/evaluate.ipynb` (evaluation)
- **Scripts**: `bash pretrain/run.sh [dataset] [nfold] [config]`, `bash finetune/run.sh <pretrained_path> [dataset] [nfold] [config]`
- **FDA inference**: `vcms/fda/` — 5-step pipeline (query drugs → embeddings → merge → ensemble predict → postprocess deg/rank)

## Config System

`src/config.py:load_config(name)` loads `<name>.yaml` from the given path (relative to `PROJECT_ROOT`) and resolves all relative paths to absolute.

**Path resolution** (`resolve_paths()`):
- `dataset.*_data_root` keys → resolved relative to `PROJECT_ROOT` (xpert/)
- `model.ATTN.ppi_gene_vector_path` → resolved relative to `PROJECT_ROOT`
- `drug_embeddings.*` keys (unimol, KPGT, morgan, smi, HG) → resolved relative to `PROJECT_ROOT`

**Output directory** (`get_output_dir()`):
- `pretrain` → `pretrain/output/<dataset>/<timestamp>/`
- `finetune` → `finetune/output/<pretrain_dataset>-<dataset>-<nfold>/<timestamp>/`
- `infer` → `metrics/output/<timestamp>/`

## Model Summary

- **~25M parameters**, hidden_size=256, 8 attention heads
- **Base encoder** (SA×4): encodes unperturbed cell state, gene-gene interactions
- **Pert encoder** (CA+SA+SA+CA): cross-attention fuses drug with cell context
- **Outputs**: x_pert, x_ctl, x_deg (= x_pert − x_ctl)
- **Loss**: 0.2·L_pert + 0.003·L_ctl + 0.2·L_deg + 1.0·L_pcc
- **Drug tokens**: UniMol molecular features (512-dim → linear proj to 256) + HG biological embedding + dose/time learnable embeddings
- **Cell tokens**: binned expression embedding + PPI gene network positional embedding
- **Attention**: Uses PyTorch 2.0+ `F.scaled_dot_product_attention` (no external flash-attn dependency)

## Key Changes from History

1. **Removed HG model training/inference** (`pretrain_hg.py`, `utils_hg.py`, `HG_data/`) — loads pre-computed HG embeddings directly
2. **Workflow-oriented directory structure** (`pretrain/`, `finetune/`, `metrics/`) replacing flat `models/`, `datasets/`, `scripts/` layout
3. **Package-based source** (`src/` with `__init__.py` lazy imports) replacing root-level module imports
4. **CLI-driven hyperparameters** — all training args overridable via CLI (`--train_lr`, `--batch_size`, etc.), not just YAML
5. **Env-var-driven shell scripts** (`pretrain/run.sh`, `finetune/run.sh`) replacing per-scenario hardcoded scripts
6. **`F.scaled_dot_product_attention`** replacing `flash_attn` library — removes external dependency
7. **`torch.amp`** replacing `torch.cuda.amp` — PyTorch 2.0+ compatible
8. **Unified `drug_embeddings/` directory** in config (`config['drug_embeddings']`) replacing scattered paths in `config['dataset']`
9. **Optimized `src/dataset.py`** with pre-extracted numpy arrays for faster `__getitem__`
