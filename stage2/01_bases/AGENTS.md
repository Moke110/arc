# AGENTS.md — stage2/01_bases

Routing for the standalone Stage 2 foundation; commands and contracts are in
[`README.md`](README.md). ADR 0061 defines the frozen gene and evaluation semantics.

- Identity: `execution/bases/identity.py`, `execution/build_identity.py` and
  `output/identity-v1/`; native evidence tables are owned within this step.
- Raw preparation: `execution/build_native.py`, `execution/bases/native.py` and
  `execution/bases/raw.py`. Read the repository execution-environment rules before
  data access. External `source/` inputs are immutable.
- Shared expression: `execution/build_shared.py` and `execution/bases/raw.py`;
  raw files are bound independently from the identity tables. Current rebuilds
  use `shared-v2` and `training-v2`; historical v1 seals stay separate.
- Training corpora: `execution/build_corpora.py`, `execution/verify_training_source.py`
  and `execution/bases/corpora.py`; bounded raw-count replay verifies all three layers.
- Evaluator-only truth: `execution/build_significance.py` and
  `execution/bases/significance.py`; ordinary training readers open only the h5ad.
- HVG/metrics: `execution/hvg_metrics.md` routes the train-only selector,
  `build_hvg.py`, `evaluate_predictions.py` and `replay_pr.py`.
