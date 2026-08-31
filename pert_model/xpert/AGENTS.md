# AGENTS.md

XPert: dual-branch Transformer (~25M params) predicting drug-induced transcriptional perturbation (978 L1000 landmark genes). Input: drug features (UniMol 3D + HG biological + dose/time tokens) and cell state (binned expression + PPI gene vectors); output: x_pert, x_ctl, x_deg = x_pert − x_ctl. Pure library — finetune/predict artifacts are owned by the workflow.

## Structure

- `src/` — the model library: `train.py` (CLI entry), `model_xpert.py` (XPertNet + AttnEncoder), `model_utils.py` (attention/embedding blocks), `dataset.py` (MyDataset), `utils.py` (dataloaders, losses, EarlyStopping), `config.py` (YAML + root-aware path resolution), `metrics.py`/`eval_metrics.py`, `predict.py` (FDA prediction), `vis_*.py`.
- `pretrain/` — the pretrain recipe (`run.sh`); config/data/weights moved to the consuming workflow per ADR 0008.
- `metadata/` — gene/drug reference (`l1000_gene_info_978.csv`, `PPI_gene_vector_128d.npy`).
- `metrics/` — evaluation scripts/notebooks. `reproducing/` — figure notebooks.
- `finetune/`, `predict/`, drug embedding files moved to the workflow per ADR 0002; pretrain weight/YAML/h5ads moved per ADR 0008.

## Key files

- `src/train.py` — arg parsing, train/validate/infer/test loops; loss 0.2·L_trt + 0.003·L_ctl + 0.2·L_deg + 1.0·L_pcc; GradScaler, early stopping on PCC (patience 50), lr scheduler. Test-mode weight load is **tolerant** (mirrors the train-mode finetune load): keys with mismatched name/shape are skipped and logged, so a pretrain checkpoint whose aux `class_fc` head drifted across library versions can still be zero-shot-evaluated (used by the workflows' `pretrain` benchmark column).
- `src/dataset.py` — AnnData → torch tensors; binned expression (n_bins=128), pre-extracted metadata numpy arrays.
- `src/utils.py` — `_load_data_and_drug_feat` (h5ad + drug feat array), `load_dataloader/test/infer`, UniMol feature gen.
- `src/config.py` — configs resolve against the owning project root (workflow-owned finetune configs resolve vs the workflow).
- `src/predict.py` — indexes the precomputed UniMol array by drug index (no SMILES regeneration), ensemble inference.

## Dependencies

- Embeddings loaded from the workflow (`workflow/xpert-fda-miri-v1/input/` + `01_drug_embed/output/`): UniMol `all_drugs_unimol_arr.npy` (8981×122×514, L1000 rows), HG `all_drugs_idx2HG.npy` (8981×256, author reference byte-identical), SMILES `all_drugs_idx2smi.npy`.
- PPI vectors: `metadata/PPI_gene_vector_128d.npy` (978×128), model-owned.
- Pretrain weight + pretrain YAML: `workflow/xpert-fda-miri-v1/04_xpert_pretrain/`; L1000 pretrain h5ads: `pert_data/l1000-xpert/data/`.
- Env `arc` (torch 2+, scanpy, numpy, scipy, pandas, scikit-learn, matplotlib, unimol_tools, pyyaml, tqdm).

## Entry points

- Pretrain: `bash pretrain/run.sh` (env vars: TRAIN_LR, BATCH_SIZE, NUM_EPOCHS, LOSS_WEIGHT, DRUG_FEAT, DRUG_EMBEDS). Direct: `python -m src.train --mode train --config ../../workflow/xpert-fda-miri-v1/04_xpert_pretrain/config/l1000 --dataset l1000_mdmt_full --device cuda:0`.
- Finetune / benchmark / predict: run via the workflow steps (`workflow/xpert-fda-miri-v1/05_xpert_finetune_bench/`, `06_xpert_finetune/`, `07_xpert_predict/`) — see its AGENTS.md.

## Notes

- Both pretrain and finetune data are log2 expression space — no z-score conversion.
- `trt_raw_data = ctl_raw_data` in FDA prediction mode; FDA prediction ensembles all `.pth` in the model dir.
- The finalized MIRI screen uses only the author-reference HG table (8981 rows) — no HG finetuned/extension rows; the HG table is `01_drug_embed/output/all_drugs_idx2HG.npy` byte-identical to `input/HG_drug_embeddings.npy`.
- Reproducibility: `set_random_seed` + env `XPERT_DETERMINISTIC=1` + `CUBLAS_WORKSPACE_CONFIG=:4096:8` for bit-exact finetune weights/predictions.
- Memory: UniMol array (float64, ~5.2 GiB) is memory-mapped at load (`mmap_mode='r'` in `utils.py` and `predict.py`) so parallel finetunes don't OOM WSL; finetune peak RSS ~2.1 GB. Row-wise float32 casting keeps values identical.

## DTOXS / CM splits

- The three MIRI workflows generate the fold columns under the **unified fold
  scheme** (ADR 0014): `warm-start` (10, drug-overlap — every test drug also
  in train; per-drug ~20% sample split), `cold-drug` (10, drug-disjoint — no
  test drug in train; sample-weighted ~20% of drugs, coverage-rescued),
  and the all-train full column (`dtoxs_full` / `cm_full`). `cold-cate` is
  dropped (no category/data-source splits). Values are `train`/`test`/`''`;
  per-fold test fraction in [0.10, 0.30]; no separate validation set
  (val = test when empty). Splits persisted in h5ad `.obs` columns.
- `cm` is an additive dataset name (union-CM finetune, `workflow/xpert-fda-miri-v2`): `_load_data_and_drug_feat` reads `DTOX_data_root` for `cm` exactly as for `dtoxs`, and the model treats `cm` like `dtoxs` for dose/time embeddings (additive branches; existing `dtoxs` behavior unchanged). Config dataset block and CLI dataset name are both `cm`; splits are `warm-start` + `cold-drug` + `cm_full`.
- Finetune defaults: 4 epochs, batch 32, init_epoch 0, finetune_lr 0.002, output to workflow `finetune/output/{dataset}_{split}/` (no timestamp subdir). All three workflows run 4 epochs (v2 aligned to 4 in ADR 0014).
