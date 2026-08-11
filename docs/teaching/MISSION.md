# Mission: HG & XPert Training Pipeline (v0 vs v1)

## Why

You work on the MIRI drug-repositioning pipeline. The two trained models — HG (HeteroGNN drug embeddings) and XPert (dual-branch Transformer DEG predictor) — are the engines behind every drug ranking you produce. You want to hold their full training path in your head: what data each stage consumes, what the loss optimizes, which artifacts come out, and how that differs between the legacy `xpert-fda-v0` flow and the refactored `xpert-fda-v1` workflow. This lets you explain, reproduce, and safely modify either model without breaking parity.

## Success looks like

- Draw from memory the 4 training flowcharts (HG pretrain, HG FDA finetune, XPert pretrain, XPert finetune) in `[input]+[model]->[output]` notation
- State the key hyperparameters and loss for each stage without looking
- Explain the v0→v1 orchestration diff: what moved, who owns artifacts, why parity matters
- Point to the exact file for any stage (e.g. "HG finetune is `drug_emb_model/hg/src/finetune_hg.py`", "XPert finetune v1 is workflow step 04")

## Constraints

- Timeboxed learning: single focused lesson now, revisitable reference docs after
- You already know the pipeline's data facts (11,104 drug indices, log2 space, 978 landmarks) — don't re-teach
- Flowcharts are the preferred notation

## Out of scope

- The scoring/ranking half of the workflow (steps 07–11) — that's a different lesson
- HG model internals beyond training (relation predictors, DTI score heads)
- Implementing or modifying the models
