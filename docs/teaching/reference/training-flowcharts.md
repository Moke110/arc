# HG & XPert Training — Flowcharts (ascii reference)

Notation: `[category]label` = node; `+` = parallel inputs; `->` = transformation
(training / data collection / inference). Companion to lesson
`lessons/0001-hg-xpert-training-flowcharts.html`.

## 1. HG pretrain (L1000, 8,981 drugs)
`drug_emb_model/hg/src/pretrain_hg.py` · launcher `run.sh`

```
[drug nodes]all_drug_node_feat_8981.npy (512d)
+ [gene nodes]all_gene_node_feat_19392.npy (512d)
+ [edges]PPI 901,260 gene-gene (STRING>700, weighted)
+ [edges]DTI 12,890 drug-gene
+ [edges]DDS 287,834 drug-drug (Tanimoto)
+ [model]HeteroGNN_untrained  (3xSAGEConv, 512->256, ReLU)
-> [train]InfoNCE contrastive · lr 5e-4 · NeighborLoader 10/5/3 ·
          batch 1 · temp 0.5 · num_neg 3 · early-stop on loss
-> [model]HG_l1000 = model/20260518_102121/global_fold_early_stop.pth (epoch 23)
-> [embs]phase1_drug_embeddings.npy (8981x256) == author reference
```

## 2a. HG FDA data collection (2,123 drugs)
`src/query_dti_fda.py` · `src/compute_dds_fda.py` · `src/infer_fda_hg.py`

```
[drugs]2,123 FDA SMILES
-> [query_dti_fda.py]ChEMBL API · pChEMBL>=6.0
        -> [edges]DTI_edge_index_fda_ext.npy (3,984 drug-gene)
+ [drugs]FDA Morgan fingerprints
-> [compute_dds_fda.py]Tanimoto>0.5106 · density-matched
        -> [edges]DDS_edge_index_fda_drugs.npy (7,659 drug-drug)
+ [drugs]FDA KPGT embeddings
-> [infer_fda_hg.py]Ridge regression KPGT->512
        -> [nodes]FDA drug node features (2,123x512)
```

## 2b. HG finetune (11,104 drugs)
`src/finetune_hg.py`

```
[model]HG_l1000
+ [graph]extended · 8,981 train + 2,123 FDA drug nodes · 19,392 gene nodes
+ [edges]train PPI/DTI/DDS + FDA DTI/DDS
-> [train]full-graph forward · InfoNCE ·
          DTI all + DDS 20K sampled + PPI 50K sampled ·
          lr 5e-5 · 50 epochs · guard: train_cosim >= 0.95 (check every 5 ep)
-> [embs]all_drugs_idx2HG.npy (11,104x256) · FDA/Train std ratio 1.000
```
v1 = workflow step 02 invokes this with `--embeddings-output` into its own dir.
v0 = ran directly inside `drug_emb_model/hg/`.

## 3. XPert pretrain (L1000)
`pert_model/XPert/source/pretrain/run.sh` · `pretrain/config/l1000.yaml` · `src/train.py`

```
[drugs]UniMol 3D feats (11,104x122x514)
+ [drugs]HG bio embeddings (11,104x256)          <- from HG
+ [drugs]dose / time / cell tokens
+ [cells]binned expression (n_bins=128) + PPI gene vectors (978x128)
+ [data]l1000_mdmt_full_336852.h5ad
+ [model]XPert_untrained · dual-branch Transformer ~25M ·
          trt CA+SA+SA+CA · ctl SA+SA+SA+SA
-> [train]loss = 0.2*L_trt + 0.003*L_ctl + 0.2*L_deg + 1.0*L_pcc ·
          lr 0.004 · batch 64 · GradScaler · early-stop on PCC (patience 50)
-> [model]XPert_l1000 = pretrain/output/pretrain_mdmt_full_50_epoch.pth
```

## 4. XPert finetune + FDA predict (DToxS)

### v0 · library-driven (`pert_model/XPert/`, pre-refactor git 47e1604^)
```
[model]XPert_l1000
+ [data]finetune/data/dtoxs.h5ad
+ [feats]drug_embeddings/ (UniMol+HG)
-> bash finetune/run.sh <ckpt> dtoxs dtoxs finetune/config/dtoxs
-> [train]4 ep · batch 32 · finetune_lr 0.002 · init_epoch 0 ·
          warm-start split -> final model
-> [model]finetune/output/dtoxs_warm/*_best.pth
-> bash predict/predict.sh --model_path finetune/output/dtoxs_warm \
        --ctl_csv cm_ctl_landmark.csv --drug_csv fda_drugs.csv
-> [output]xpert_fda-cm_deg.csv
```

### v1 · workflow step-wise (`stage1/xpert-fda-miri-v1/`)
```
[model]XPert_l1000
+ [data]input/dtoxs.h5ad
+ [feats]input/* + 02_drug_embeddings/output/*      (via config/dtoxs.yaml)
-> 04_finetune_xpert.py   (dtoxs all-train · 4 ep · batch 32 · lr 0.002)
-> [model]04_finetune_xpert/output/dtoxs_dtoxs/dtoxs_dtoxs_fold-1_best.pth
-> 05_benchmark_xpert.py  (warm-start 10f · cold-drug 10f · cold-cate 6f)
-> [eval]05/output/benchmark.csv + benchmark.png
-> 06_perturbed_transcriptome.py
-> [output]06/output/xpert_fda-cm_deg.csv
```

## Cheat sheet

| Stage | Model | Objective | Key params | Artifact |
|-------|-------|-----------|------------|----------|
| HG pretrain | HeteroGNN 3xSAGEConv 512->256 | InfoNCE (temp 0.5) | lr 5e-4, cap 300 ep, pat 30, num_neg 3, NB 10/5/3 | `global_fold_early_stop.pth` · 8981x256 |
| HG finetune | same HeteroGNN | InfoNCE (DTI all + DDS/PPI sampled) | lr 5e-5, 50 ep, cosim>=0.95 | `all_drugs_idx2HG.npy` 11104x256 |
| XPert pretrain | dual-branch Transformer ~25M | 0.2·L_trt+0.003·L_ctl+0.2·L_deg+1.0·L_pcc | lr 0.004, batch 64, pat 50 (PCC) | `pretrain_mdmt_full_50_epoch.pth` |
| XPert finetune | XPert (from pretrained) | same 4-term loss | 4 ep, batch 32, lr 0.002, init 0 | `dtoxs_dtoxs_fold-1_best.pth` (v1) |

## Key v0 vs v1 differences
- v0: scripts/notebooks inside `pert_model/XPert/`; library owns artifacts; final
  model from warm-start split.
- v1: XPert is a pure library; workflow steps 02/04/05/06 own weights + CSVs;
  final model from all-train `dtoxs` split; benchmark folds in step 05.
- v1 steps are idempotent (skip on cached artifacts); retraining is
  GPU-fp-nondeterministic and can break parity with archived CSVs.
