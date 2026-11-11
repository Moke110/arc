# AGENTS.md — Stage 2 shared retrieval evaluation

- [README.md](README.md) and [ADR 0059](../../docs/adr/0059-stage2-significance-retrieval.md) define the tenv-only metric contract and evaluator/predictor boundary.
- `metrics.py` owns exact group AP and bounded-memory pooled PR; every step uses this implementation.
- `build_labels.py` owns checkpointed raw-count DESeq2 label building. Dataset jobs use the repository resource limits and dedicated scheduling owner.
- `inventory.json` binds frozen predictor artifacts; `inventory.py` verifies these bindings.
- `evaluate.py` scores and replays; `publish.py` independently scores in two processes before replacing a generated release.
- `tests/` covers the scientific definitions, immutable-input guards and all five reporting entrypoints.
