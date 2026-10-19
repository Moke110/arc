"""
Model-adapter registry for the Stage 2 v1 benchmark.

Entrants are added here as their adapters land (tickets 07-11). Each entry maps
a model name to an adapter class implementing the `predict(train, test)`
contract in `_lib/adapter.py`; adapters drive their own pinned conda environment
through the runner harness (`_lib/model_runner.py`).

Until a model adapter exists it is absent here, so `04_evaluate` simply does not
score it; the checkpoint audit (`docs/entrant_audit.md`) records why.
"""

import importlib.util
from pathlib import Path

WORKFLOW_ROOT = Path(__file__).resolve().parent.parent


def _load_state():
    path = WORKFLOW_ROOT / "03_state" / "adapter.py"
    spec = importlib.util.spec_from_file_location("cm_bench_state_adapter", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.StateAdapter


MODEL_ADAPTERS = {}
try:
    MODEL_ADAPTERS["state"] = _load_state()
except Exception as exc:  # adapter import must never break baselines
    print(f"[registry] state adapter unavailable: {exc}")
