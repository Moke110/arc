# AGENTS.md

XPert: dual-branch Transformer (~25M params) predicting drug-induced transcriptional perturbation (978 L1000 landmark genes). Input: drug features (UniMol 3D + HG biological + dose/time tokens) and cell state (binned expression + PPI gene vectors); output: x_pert, x_ctl, x_deg = x_pert − x_ctl. Pure library — finetune/predict artifacts are owned by the workflow.

## Structure

- `src/` — the model library: `train.py` (CLI entry), `model_xpert.py` (XPertNet + AttnEncoder), `model_utils.py` (attention/embedding blocks), `dataset.py` (MyDataset), `utils.py` (dataloaders, losses, EarlyStopping), `config.py` (YAML + root-aware path resolution), `metrics.py`/`eval_metrics.py`, `predict.py` (FDA prediction), `vis_*.py`.
- `pretrain/` — L1000 pretraining (configs, data, `run.sh`, output weights).
- `metadata/` — gene/drug reference (`l1000_gene_info_978.csv`, `PPI_gene_vector_128d.npy`).
- `metrics/` — evaluation scripts/notebooks. `reproducing/` — figure notebooks.
- `finetune/`, `predict/`, drug embedding files moved to the workflow per ADR 0002.

## Key files

- `src/train.py` — arg parsing, train/validate/infer/test loops; loss 0.2·L_trt + 0.003·L_ctl + 0.2·L_deg + 1.0·L_pcc; GradScaler, early stopping on PCC (patience 50), lr scheduler.
- `src/dataset.py` — AnnData → torch tensors; binned expression (n_bins=128), pre-extracted metadata numpy arrays.
- `src/utils.py` — `_load_data_and_drug_feat` (h5ad + drug feat array), `load_dataloader/test/infer`, UniMol feature gen.
- `src/config.py` — configs resolve against the owning project root (workflow-owned finetune configs resolve vs the workflow).
- `src/predict.py` — indexes the precomputed UniMol array by drug index (no SMILES regeneration), ensemble inference.

## Dependencies

- Embeddings loaded from the workflow (`miri_reposition/xpert-fda-v1/input/` + `02_drug_embeddings/output/`): UniMol `all_drugs_unimol_arr.npy` (11104×122×514), HG `all_drugs_idx2HG.npy` (11104×256), SMILES `all_drugs_idx2smi.npy`.
- PPI vectors: `metadata/PPI_gene_vector_128d.npy` (978×128), model-owned.
- Env `arc` (torch 2+, scanpy, numpy, scipy, pandas, scikit-learn, matplotlib, unimol_tools, pyyaml, tqdm).

## Entry points

- Pretrain: `bash pretrain/run.sh` (env vars: TRAIN_LR, BATCH_SIZE, NUM_EPOCHS, LOSS_WEIGHT, DRUG_FEAT, DRUG_EMBEDS). Direct: `python -m src.train --mode train --config pretrain/config/l1000 --dataset l1000_mdmt_full --device cuda:0`.
- Finetune / benchmark / FDA predict: run via the workflow steps (`miri_reposition/xpert-fda-v1/04_finetune_xpert/`, `05_benchmark_xpert/`, `06_perturbed_transcriptome/`) — see its AGENTS.md.

## Notes

- Both pretrain and finetune data are log2 expression space — no z-score conversion.
- `trt_raw_data = ctl_raw_data` in FDA prediction mode; FDA prediction ensembles all `.pth` in the model dir.
- HG covers 11,104 drugs: training rows 0–8980 = author reference; FDA rows 8981–11103 = fine-tuned + quantile norm mapping (std ratio 1.000).
- Reproducibility: `set_random_seed` + env `XPERT_DETERMINISTIC=1` + `CUBLAS_WORKSPACE_CONFIG=:4096:8` for bit-exact finetune weights/predictions.
- Memory: UniMol array (float64, ~5.2 GiB) is memory-mapped at load (`mmap_mode='r'` in `utils.py` and `predict.py`) so parallel finetunes don't OOM WSL; finetune peak RSS ~2.1 GB. Row-wise float32 casting keeps values identical.

## DTOXS splits

- cold-drug: 10-fold leave-one-drug-out; cold-cate: 6-fold leave-one-category-out; warm-start: 10-fold drug-stratified random (evaluation-only); dtoxs: 1-fold all-train (final model). Splits persisted in h5ad `.obs` columns; no separate validation set (val = test when empty).
- Finetune defaults: 4 epochs, batch 32, init_epoch 0, finetune_lr 0.002, output to workflow `finetune/output/dtoxs_{split}/` (no timestamp subdir).
