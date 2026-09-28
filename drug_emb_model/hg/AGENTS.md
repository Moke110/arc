# AGENTS.md

File-content index for `drug_emb_model/hg/` — HeteroGNN drug biological-embedding library (dormant; no canonical workflow consumes it).

- `src/` — model + training/inference code (see below).
- `run.sh` — training launcher.
- `run.ipynb` — interactive notebook.
- `src/model_hg.py` — HeteroGNN model definition.
- `src/pretrain_hg.py` — contrastive pretraining entry.
- `src/finetune_hg.py` — extended-graph fine-tune.
- `src/infer_fda_hg.py` — FDA-extension inference.
- `src/query_dti_fda.py` — ChEMBL DTI edge collection.
- `src/compute_dds_fda.py` — Tanimoto DDS edge computation.
- `src/gene_id_mapping.py` — gene ID mapping helpers.
- `src/merge_dti_batches.py` — merge DTI batch outputs.
- `src/utils_hg.py` — contrastive loss, sampling, embedding utils.
- `src/__init__.py` — package init.

Out-of-dir: author HG reference in the workflow `input/`; archived FDA data/weights in `archived/hg-fda-extension-archive/`.
