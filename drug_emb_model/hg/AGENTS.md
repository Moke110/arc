# AGENTS.md

Heterogeneous Graph Neural Network (HeteroGNN) for drug biological embedding generation. Trains on a drug-gene heterogeneous graph via contrastive learning (InfoNCE) to produce 256-dim biological drug embeddings consumed by XPert.

**Pure static library (ADR 0008)**: this directory holds only `src/` code, `run.sh`, `run.ipynb`, and this doc. All data/weights moved out — the author HG reference now lives in the workflow `input/`, and the orphaned FDA-extension data/weights are archived (gitignored) at `reference/hg-fda-extension-archive/` (layout mirrors the origin component paths). The model is dormant: no canonical workflow consumes it.

## Structure

- `src/` — `model_hg.py` (HeteroGNN: 3 HeteroConv layers with SAGEConv, relation MLP heads), `pretrain_hg.py` (training entry), `finetune_hg.py` (extended-graph fine-tune, low LR, norm mapping), `infer_fda_hg.py`, `query_dti_fda.py` (ChEMBL DTI), `compute_dds_fda.py` (Morgan Tanimoto DDS edges), `gene_id_mapping.py`, `merge_dti_batches.py`, `utils_hg.py` (ContrastiveLoss, negative sampling, get_embeddings).
- `run.sh` — training launcher; `run.ipynb` — interactive notebook.

## Preserved embeddings

| File | Shape | Status |
|------|-------|--------|
| `workflow/xpert-fda-miri-v1/input/HG_drug_embeddings.npy` | (8981, 256) | Author L1000 reference, workflow-owned |
| `workflow/xpert-fda-miri-v1/01_drug_embed/output/all_drugs_idx2HG.npy` | (8981, 256) | Workflow step-01 copy, byte-identical to the author reference |

The FDA-extension rows (8981–11103) and the HG checkpoints/embeddings that produced them are archived at `reference/hg-fda-extension-archive/` (data/, embeddings/, model/, figshare/).

## Key workflows

```bash
# Training from scratch (8,981 drugs) — needs the archived data/ restored:
#   cp -r reference/hg-fda-extension-archive/data .
bash run.sh                                # defaults: 3 layers, hidden 256, lr 5e-4, epochs 300, patience 30, temperature 0.5, num_neg 3
```

## Model & data (archived)

- Model: 3 HeteroConv layers (SAGEConv per edge type), gene/drug inputs 512→256, ReLU after each layer, 3 relation predictors (Linear(512→1)) for DTI/DDS/PPI; ~3M params.
- Graph: PPI 901,260 weighted gene-gene (STRING >700); DTI 12,890 + 3,984 FDA (ChEMBL pChEMBL≥6.0); DDS 287,834 + 7,659 FDA (Tanimoto>0.5106); drug nodes 8,981 + 2,123 = 11,104; gene nodes 19,392 (978 L1000 landmarks).
- Loss: InfoNCE contrastive (temperature 0.5); positives = real edges, negatives = sampled non-edges (3–5 per positive).

## Status

- Pretraining complete (2026-05-18, epoch 23, early-stopped); FDA data collection complete (3,984 DTI + 7,659 DDS edges); fine-tuning complete (2026-06-10, 50 epochs). All artifacts preserved in the archive for provenance.
