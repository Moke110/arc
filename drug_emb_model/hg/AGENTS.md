# AGENTS.md

Heterogeneous Graph Neural Network (HeteroGNN) for drug biological embedding generation. Trains on a drug-gene heterogeneous graph via contrastive learning (InfoNCE) to produce 256-dim biological drug embeddings consumed by XPert.

## Structure

- `src/` — `model_hg.py` (HeteroGNN: 3 HeteroConv layers with SAGEConv, relation MLP heads), `pretrain_hg.py` (training entry), `finetune_hg.py` (extended-graph fine-tune, low LR, norm mapping), `infer_fda_hg.py`, `query_dti_fda.py` (ChEMBL DTI), `compute_dds_fda.py` (Morgan Tanimoto DDS edges), `gene_id_mapping.py`, `merge_dti_batches.py`, `utils_hg.py` (ContrastiveLoss, negative sampling, get_embeddings).
- `data/` — precomputed node features and edge indices (PPI/DTI/DDS, gene+FDA drug nodes, 978 landmark gene indices, FDA drug list).
- `model/` — pretrain + finetune checkpoints and phase-1 drug embeddings. `embeddings/article/HG_drug_embeddings.npy` — author reference (8981×256).
- `run.sh` — training launcher; `run.ipynb` — interactive notebook.

## Preserved embeddings

| File | Shape | Status |
|------|-------|--------|
| `miri_reposition/xpert-fda-v1/02_drug_embeddings/output/all_drugs_idx2HG.npy` | (11104, 256) | Active, workflow-generated |
| `embeddings/article/HG_drug_embeddings.npy` | (8981, 256) | Author L1000 reference |

Training rows (0–8980) identical to the author reference; FDA rows (8981–11103) from fine-tuning on the extended graph + quantile-interpolation norm mapping (FDA/Train std ratio 1.000, DTI signal 2.48x vs random). Per ADR 0002, HG does not write into a workflow dir at runtime: `finetune_hg.py --embeddings-output` writes the all-drug table, and workflow step 02 invokes it with the workflow path.

## Key workflows

```bash
# Training from scratch (8,981 drugs)
bash run.sh                                # defaults: 3 layers, hidden 256, lr 5e-4, epochs 300, patience 30, temperature 0.5, num_neg 3

# FDA drug embedding generation (full pipeline)
python -m src.query_dti_fda --start 0 --end 2123 --min-pchembl 6.0
python -m src.compute_dds_fda
python -m src.finetune_hg --device cuda:0 --lr 5e-5 --epochs 50 \
    --update-xpert-embeddings --embeddings-output embeddings/all_drugs_idx2HG.npy

# Resume
bash run.sh --global_saved_model model/<timestamp>/global_fold_early_stop.pth
```

## Model & data

- Model: 3 HeteroConv layers (SAGEConv per edge type), gene/drug inputs 512→256, ReLU after each layer, 3 relation predictors (Linear(512→1)) for DTI/DDS/PPI; ~3M params.
- Graph: PPI 901,260 weighted gene-gene (STRING >700); DTI 12,890 + 3,984 FDA (ChEMBL pChEMBL≥6.0); DDS 287,834 + 7,659 FDA (Tanimoto>0.5106); drug nodes 8,981 + 2,123 = 11,104; gene nodes 19,392 (978 L1000 landmarks).
- Loss: InfoNCE contrastive (temperature 0.5); positives = real edges, negatives = sampled non-edges (3–5 per positive).
- Fine-tune: full-graph training, low LR (5e-5) preserves training embeddings (cosim>0.95), then quantile-interpolation norm mapping.

## Status

- Pretraining complete (2026-05-18, epoch 23, early-stopped); FDA data collection complete (3,984 DTI + 7,659 DDS edges); fine-tuning complete (2026-06-10, 50 epochs). Checkpoints preserved for both.
