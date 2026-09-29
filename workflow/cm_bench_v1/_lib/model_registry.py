"""
Model-adapter registry for the Stage 2 v1 benchmark.

Entrants are added here as their adapters land (tickets 07-11). Each entry maps
a model name to an adapter class implementing the `predict(train, test)`
contract in `_lib/adapter.py`; adapters drive their own pinned conda environment
through the runner harness (`_lib/model_runner.py`).

Until a model adapter exists it is absent here, so `04_evaluate` simply does not
score it; the checkpoint audit (`docs/entrant_audit.md`) records why.
"""

MODEL_ADAPTERS = {}
