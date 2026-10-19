# HG & XPert Training Resources

High-trust sources for the training pipeline. Primary sources are in-repo (this is your own codebase); the two literature items ground the models.

## Knowledge

### HG (HeteroGNN)

- [AGENTS.md: `drug_emb_model/hg/AGENTS.md`](../../drug_emb_model/hg/AGENTS.md)
  The fastest overview of the HG component. Use for: data graph facts (edge counts, node counts), preserved-embeddings table, key workflow commands, status.
- [`src/pretrain_hg.py`](../../drug_emb_model/hg/src/pretrain_hg.py)
  Pretrain entrypoint. Use for: graph construction (edge types, HeteroData), NeighborLoader setup, InfoNCE contrastive training loop, early stopping.
- [`src/finetune_hg.py`](../../drug_emb_model/hg/src/finetune_hg.py)
  FDA-extension finetune. Use for: full-graph forward, low-LR design, DTI/DDS/PPI per-edge loss, cosim≥0.95 guard, norm-mapping rationale, CLI (`--update-xpert-embeddings --embeddings-output`).
- [`src/infer_fda_hg.py`](../../drug_emb_model/hg/src/infer_fda_hg.py)
  FDA node feature generation (Ridge KPGT→512) and extended-graph build. Use for: how FDA drugs enter the graph.
- [`src/query_dti_fda.py`](../../drug_emb_model/hg/src/query_dti_fda.py) + [`src/compute_dds_fda.py`](../../drug_emb_model/hg/src/compute_dds_fda.py)
  FDA edge collection: ChEMBL API DTI (pChEMBL≥6.0) and Morgan/Tanimoto DDS (density-calibrated).
- [`run.sh`](../../drug_emb_model/hg/run.sh)
  Canonical pretrain invocation with exact defaults (3 layers, 256 hidden, lr 5e-4, 300 epochs, patience 30, temp 0.5, num_neg 3).

### XPert

- [AGENTS.md: `pert_model/xpert/AGENTS.md`](../../pert_model/xpert/AGENTS.md)
  Overview of XPert as a pure library. Use for: model structure, loss weights, DTOXS splits, entry points, reproducibility notes.
- [`pretrain/run.sh`](../../pert_model/xpert/pretrain/run.sh)
  Pretrain launcher with defaults (lr 0.004, batch 64, 2500 epochs cap, loss 0.2,0.003,0.2,1, patience 50, drug embeds `unimol,HG`).
- [`stage1/xpert-fda-miri-v1/04_xpert_pretrain/config/l1000.yaml`](../../stage1/xpert-fda-miri-v1/04_xpert_pretrain/config/l1000.yaml)
  Pretrain recipe (workflow-owned); data roots point into `pert_data/l1000-xpert/data/` + the workflow `input/`.
- [`src/train.py`](../../pert_model/xpert/src/train.py)
  Train/finetune/test/infer loops. Use for: exact loss composition, GradScaler, early stopping on PCC, finetune param handling.
- [`src/model_xpert.py`](../../pert_model/xpert/src/model_xpert.py)
  Dual-branch Transformer (trt: CA+SA+SA+CA, ctl: SA+SA+SA+SA). Use for: how drug + cell inputs fuse.

### Workflows (v0 vs v1)

- [AGENTS.md: `stage1/xpert-fda-miri-v1/`](../../stage1/xpert-fda-miri-v1/AGENTS.md)
  Use for: workflow layout, entrypoints, data facts/gotchas, ADR 0002/0008 (artifacts owned by workflow, components pure libraries).
- [`stage1/xpert-fda-miri-v1/01_drug_embed/drug_embeddings.py`](../../stage1/xpert-fda-miri-v1/01_drug_embed/drug_embeddings.py)
  Step 01: copies the workflow-owned author HG reference (`input/HG_drug_embeddings.npy`) with a byte-parity gate.
- [`stage1/xpert-fda-miri-v1/06_xpert_finetune/xpert_finetune.py`](../../stage1/xpert-fda-miri-v1/06_xpert_finetune/xpert_finetune.py) + [`config/dtoxs.yaml`](../../stage1/xpert-fda-miri-v1/06_xpert_finetune/config/dtoxs.yaml)
  Step 06: all-train dtoxs finetune, artifacts moved into step output.
- [`stage1/xpert-fda-miri-v1/05_xpert_finetune_bench/xpert_finetune_bench.py`](../../stage1/xpert-fda-miri-v1/05_xpert_finetune_bench/xpert_finetune_bench.py)
  Step 05: warm-start / cold-drug / cold-cate benchmark folds.
- [`stage1/xpert-fda-miri-v1/07_xpert_predict/perturb_transcriptome.py`](../../stage1/xpert-fda-miri-v1/07_xpert_predict/perturb_transcriptome.py)
  Step 07: FDA prediction → `xpert_fda-cm_deg.csv`.
- v0 finetune/predict (pre-refactor, viewable via git):
  `git show 47e1604^:pert_model/xpert/finetune/run.sh` and `git show 47e1604^:pert_model/xpert/predict/predict.sh`
  The old library-owned orchestration that v1 replaced.
- [`stage1/xpert-fda-miri-v0/scripts/run_all.py`](../../stage1/xpert-fda-miri-v0/scripts/run_all.py)
  The scoring half that consumed the v0 DEG CSV (WCS+RGES / Zhang / pathway RGES).

### Literature

- [`reference/markdowns/2026_Yang_XPert.md`](../../reference/markdowns/2026_Yang_XPert.md) (PDF: [`reference/PDFs/2026_Yang_XPert.pdf`](../../reference/PDFs/2026_Yang_XPert.pdf))
  The XPert paper (transcript to md). Grounds the dual-branch design and L1000 prediction task.
- [`reference/DToxS_article.md`](../../reference/DToxS_article.md)
  DToxS (cardiotoxicity hPSC-CM) dataset. Grounds what the finetune task measures.
- [`reference/PDFs/2024_Zheng_TranSiGen.pdf`](../../reference/PDFs/2024_Zheng_TranSiGen.pdf) — the HG/contrastive-graph lineage.
- [ADR 0002: components-as-libraries](../../docs/adr/0002-components-as-libraries.md)
  Why v1 moved finetune/predict artifacts out of the libraries and into the workflow.

## Wisdom (Communities)

- None yet. This is a solo research codebase; skip community delegation unless a question needs an ML peer review (then: r/bioinformatics, r/MLQuestions).
