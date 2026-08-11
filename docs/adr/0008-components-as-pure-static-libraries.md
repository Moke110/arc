# Components as pure static libraries (extends ADR 0002)

Components (`pert_model/xpert/`, `drug_emb_model/hg/`) ship library code (`src/`) and small static metadata only. All changeable data, trained weights, and dataset h5ads live either in the consuming workflow or in the appropriate dataset home (`pert_data/`). A component does not write into a workflow's directory at runtime.

## Context

ADR 0002 made workflows own the artifacts they consume, but the components still shipped training data and weights (`pert_model/xpert/pretrain/` held the L1000 h5ads and pretrain checkpoints; `drug_emb_model/hg/` held its own checkpoints, phase-1 embeddings, FDA-extension graph caches, and figshare zips). With the finalized MIRI workflow reading the author HG table and pretrained XPert weight directly out of component directories, neither the component nor the workflow was self-contained.

## Decision

- **XPert stays code + metadata**: `src/`, `metadata/` (gene info, PPI vectors), `metrics/`, `reproducing/`, and the pretrain `run.sh` recipe stay in the library. The used pretrain checkpoint and the pretrain YAML move into the workflow (`workflow/xpert-fda-miri-v1/03_xpert_pretrain/`); the two unused pretrain checkpoints are deleted; the three L1000 pretrain h5ads move to the dataset home `pert_data/l1000/`.
- **HG stays code only**: `src/`, `run.sh`, `run.ipynb`, `AGENTS.md` remain. The author HG reference moves into the workflow `input/`; the orphaned FDA-extension data/weights (graph caches, KPGT embeddings, HG checkpoints, figshare zips) move to the gitignored archive `reference/hg-fda-extension-archive/` (paths mirror the origin component layout for provenance).
- **The finalized workflow owns all its inputs**: `input/` is a real directory holding the frozen UniMol array, idx→SMILES mapping, DToxS h5ad, FDA drug list, curated MIRI pathways, and the author HG reference; the pretrain weight + both YAML recipes (pretrain + finetune) live in the workflow.
- `pert_data/dtoxs/` and `disease_data/miri/{GSE}/` remain shared static inputs.

## Consequences

- A component is a stable, pure library; the workflow is self-contained on drug data and training configuration.
- HG is a dormant library — its `run.sh` still points at the archived `data/`; only `src/` code is authoritative.
- The finalized workflow never reads an embedding/weight asset from a component directory (verified by the library-purity and path invariants).
