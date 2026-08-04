# CLAUDE.md

## Purpose

XPert is a dual-branch transformer (~25M params) that predicts drug-induced transcriptional perturbation responses (978 L1000 landmark genes). It takes drug molecular features (UniMol + HG biological embeddings + dose/time tokens) and cell state (binned expression + PPI gene vectors), and outputs predicted perturbed expression (x_pert), control expression (x_ctl), and expression change (x_deg = x_pert - x_ctl).

## Structure

```
xpert/
├── src/                       # All source code (the model library)
│   ├── train.py               # CLI entry point: arg parsing, train/validate/infer/test loops, main()
│   ├── model_xpert.py         # XPertNet (main model) + AttnEncoder (dual-branch transformer stack)
│   ├── model_utils.py         # Attention (SA, CA), Embeddings (cell, drug, unimol), Encoder/Decoder blocks, FFN
│   ├── dataset.py             # MyDataset: AnnData -> PyTorch tensors (lazy __getitem__, pre-extracted numpy arrays)
│   ├── utils.py               # Data loading (load_dataloader), loss functions (MSE, PCC, cosine), EarlyStopping, UniMol feat gen
│   ├── config.py              # YAML config loader, root-aware path resolver, output directory helper
│   ├── metrics.py             # Core metrics: MSE/PCC/Precision@K + get_metrics_new + get_metrics_infer (raw + deg variants)
│   ├── eval_metrics.py        # Extended metrics: R2, MMD, Wasserstein
│   ├── predict.py             # FDA drug perturbation prediction: indexes precomputed UniMol array
│   ├── vis_gene_attn.py       # Gene attention heatmaps/boxplots/scatter
│   ├── vis_molecule.py        # Molecule attention visualization with RDKit
│   └── __init__.py            # Lazy imports for src package
├── pretrain/                  # Step 1: L1000 pretraining (model-owned)
│   ├── config/l1000.yaml      # Pretrain config: batch 64, lr 0.004, hidden 256, CA+SA+SA+CA / SA+SA+SA+SA
│   ├── data/                  # Input h5ad files
│   ├── output/                # Pretrained model weights
│   └── run.sh                 # Launch script (env-var-driven: DATASET, SPLIT_METHOD, CONFIG, etc.)
├── metadata/                  # Gene/drug reference info
│   ├── l1000_gene_info_978.csv # L1000 landmark gene symbols, entrez IDs, gene_idx
│   └── PPI_gene_vector_128d.npy # PPI gene network positional embeddings (978 x 128)
├── metrics/                   # Step 4: Evaluation
│   ├── evaluate.ipynb         # Interactive evaluation notebook
│   ├── evaluate_dtoxs.py      # DTOXS evaluation: 4 model groups, 14 deg-only metrics
│   ├── old-hg/                # Evaluation results with previous HG embeddings
│   └── new-hg/                # Evaluation results with new HG embeddings, 4-epoch finetune
└── reproducing/               # Figure reproduction notebooks
```

Note: `finetune/`, `predict/`, and the drug embedding files moved to the workflow
(`miri_reposition/xpert-fda-v1/`) per ADR 0002 — components are libraries and the
workflow owns artifact generation. The pretrain config references the workflow's
embedding files and the finetune/predict stages live in the workflow.

## Key Files

- src/model_xpert.py -- XPertNet (forward: 10-tuple -> x_pert, x_ctl, x_deg), AttnEncoder (CA+SA layer stack). Drug branch: UniMol feat + HG embed + dose/time tokens. Cell branch: binned expr embed + PPI gene vector.
- src/model_utils.py -- Building blocks: SelfAttention (F.scaled_dot_product_attention), CrossAttention, cell_Embeddings (expr bins + PPI pretrained), unimol_Embeddings (512-dim atom -> 256-dim), Encoder, crossEncoder, FeedForward, LayerNorm, SublayerConnection.
- src/train.py -- arg_parse() for all CLI args, main() orchestration, train()/validate()/infer()/test() loops. Loss: 0.2*L_trt + 0.003*L_ctl + 0.2*L_deg + 1.0*L_pcc. Supports GradScaler, early stopping (PCC), lr_scheduler, resume/pretrained loading.
- src/dataset.py -- MyDataset(Dataset): loads AnnData, digitizes expression into n_bins, pre-extracts all metadata columns as numpy arrays. __getitem__ returns 10-tuple (or 14-tuple for infer mode).
- src/utils.py -- _load_data_and_drug_feat (AnnData h5ad + drug feat array/dict), load_dataloader/test/infer, loss functions (mse_loss_ls_sum, pcc_loss_sum, cos_loss_sum), EarlyStopping, UniMol feature generation from SMILES.
- src/config.py -- load_config(name) loads YAML; config paths resolve against the owning project root (so a workflow-owned finetune config resolves relative to the workflow); get_config_root(), get_output_dir().
- src/predict.py -- Standalone FDA drug prediction: indexes the precomputed UniMol array by drug index (no SMILES regeneration), ensemble inference, no DataLoader dependency.
- src/metrics.py -- MSE, RMSE, MAE, Pearson, Spearman, Precision@K (positive/negative). get_metrics_new (raw + deg), get_metrics_infer (raw + deg with per-sample values).
- predict/predict.sh -- Bash wrapper for FDA prediction (now in the workflow: miri_reposition/xpert-fda-v1/predict/predict.sh).
- predict/build_baseline.py -- Reads DTOXS raw data, computes consensus baseline across 206 control wells (now in the workflow).
- finetune/create_dtoxs_splits.py -- Generates dtoxs_full obs column in dtoxs.h5ad (now in the workflow).
- metrics/evaluate_dtoxs.py -- DTOXS evaluation: 3 model groups (L1000, cold-drug, cold-cate) across 14 deg-only metrics, bar charts with error bars.

## Dependencies

- Pre-computed embeddings (loaded from the workflow's `miri_reposition/xpert-fda-v1/data/embeddings/`, referenced by the pretrain config via relative paths):
  - UniMol: all_drugs_unimol_arr.npy (11104 x 122 x 514) -- generated by unimol_tools UniMolRepr
  - HG: all_drugs_idx2HG.npy (11104 x 256) -- final full-drug embeddings, training rows = author reference, FDA rows = HG fine-tuning + quantile norm mapping (see drug_emb_model/hg/)
  - SMILES mapping: all_drugs_idx2smi.npy (11104 dict)
- PPI vectors: metadata/PPI_gene_vector_128d.npy (978 x 128) -- retained in xpert (model-owned)
- Python packages: torch (2.0+), scanpy, numpy, scipy, pandas, scikit-learn, matplotlib, unimol_tools, pyyaml, tqdm
- Conda env: arc (Python 3.13), from root environment.yml
- First 8981 drugs (indices 0-8980) seen during L1000 pretraining. Indices 8981-11103 are FDA-extended drugs with embeddings but not seen in training.

## Entry Points

- Pretrain: bash pretrain/run.sh [dataset] [split_method] [config]
  Defaults: l1000_mdmt_full, pretrain_split, pretrain/config/l1000
  Env vars: TRAIN_LR, BATCH_SIZE, NUM_EPOCHS, LOSS_WEIGHT, DRUG_FEAT, DRUG_EMBEDS

- Direct CLI: python -m src.train --mode train --config pretrain/config/l1000 --dataset l1000_mdmt_full --device cuda:0

- Finetune / FDA prediction / split creation: now live in the workflow
  (`miri_reposition/xpert-fda-v1/finetune/`, `predict/`); run via the numbered
  stages or `screen.ipynb`. See miri_reposition/xpert-fda-v1/CLAUDE.md.

- Evaluation: python metrics/evaluate_dtoxs.py or metrics/evaluate.ipynb

## Model Summary

- ~25M parameters, hidden_size=256, 8 attention heads per layer
- Drug embedding: UniMol atom features (512-dim) projected to 256-dim + HG biological embedding (256-dim) replacing first token + dose/time learnable embeddings (for mdmt datasets)
- Cell embedding: binned expression embedding (n_bins=128 vocab) + PPI gene network positional embedding (128-dim or linear-projected)
- Pert branch (trt_structure): CA + SA + SA + CA -- cross-attention fuses drug info into cell context
- Ctl branch (ctl_structure): SA + SA + SA + SA -- self-attention over cell state only
- Output heads: three separate MLPs (256 -> 64 -> 1) producing x_pert, x_ctl, x_deg = (x_pert - x_ctl)
- Loss: weighted sum of MSE(x_pert), MSE(x_ctl), MSE(x_deg), and (1 - PCC(x_deg))
  Default weights: [0.2, 0.003, 0.2, 1.0]
- Attention: F.scaled_dot_product_attention (PyTorch 2.0+) with dropout; traditional softmax path for output_attention mode
- Training: Adam, lr_scheduler (1.0 for first 40 epochs, 0.5 after), early stopping on val PCC loss (patience=50), warmup period (init_epoch=70) with batch-normalized loss before full weighted loss backprop

## Config System

- YAML files in pretrain/config/ (model-owned) and finetune/config/ (workflow-owned, in miri_reposition/xpert-fda-v1/finetune/)
- load_config(name) loads {name}.yaml; config paths resolve against the config's owning root (workflow finetune config resolves relative to the workflow)
- resolve_paths() makes all relative paths absolute: dataset data roots, ppi_gene_vector_path, drug_embeddings.*
- Output dirs: pretrain/output/<dataset>_<split>/<timestamp>/ (model-owned); finetune outputs follow the config root (workflow's finetune/output/)
- All training hyperparameters overridable via CLI args (--train_lr, --batch_size, etc.)
- Split methods defined in config YAML: pretrain_split, cold-drug, cold-cate

## DTOXS Splits

- cold-drug: 10-fold leave-one-drug-out (one drug per category held out per fold)
- cold-cate: 6-fold leave-one-category-out
- warm-start: 10-fold drug-stratified random, no validation set, evaluation-only. Each sample in >=1 test fold.
- dtoxs: 1-fold all-train (all samples used for training, no validation/test). For final prediction model.
- Splits persisted in h5ad .obs columns (cold-drug_fold-1..10, cold-cate_fold-1..6, warm-start_fold-1..10, dtoxs_full)
- No separate validation set in cold-drug/cold-cate/warm-start/dtoxs: val_data = test_data when val is empty
- Only deg-variant metrics retained (MAE_deg, MSE_deg, Pearson_deg, Spearman_deg, R2_deg, Precision@K_deg)
- Finetuning: 4 epochs default (NUM_EPOCHS=4 in run.sh), pretrained_model=pretrain_mdmt_full_50_epoch.pth. Output saved directly to the workflow's finetune/output/dtoxs_{split}/ (no timestamp subdirectory).

## Notes

- Both pretrain (L1000) and finetune (DTOXS) data are in log2 expression space -- no z-score conversion.
- trt_raw_data = ctl_raw_data in FDA prediction mode (actual treatment expression unknown, only baseline available).
- FDA prediction uses ensemble averaging across all .pth files in model directory.
- HG embeddings cover 11,104 drugs. Training rows (0-8980) = author reference. FDA rows (8981-11103) = HG fine-tuning + quantile norm mapping. FDA/Train std ratio = 1.000.
- Drug embedding features: all arrays use same index scheme (idx->smi from all_drugs_idx2smi.npy), consumed by drug_idx from the workflow's data/embeddings/.
- predict.py indexes the precomputed UniMol array by drug index (no SMILES regeneration), so finetune and predict consume identical drug features.
- Model evaluation: evaluate_dtoxs.py runs inference on 4 model groups (L1000 pretrain, warm-start finetune, cold-drug finetune, cold-cate finetune) and reports 14 deg-only metrics (MSE, RMSE, MAE, Pearson, Spearman, R2, Precision@10/20/50/100 Positive/Negative). Results saved as CSV + PNG bar charts to metrics/dtoxs/.
