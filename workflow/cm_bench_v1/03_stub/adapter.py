"""
Stub adapter: proof that the runner harness works end-to-end on CPU.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from model_runner import SubprocessAdapter  # noqa: E402

ENTRY = Path(__file__).resolve().parent / "entry.py"


class StubAdapter(SubprocessAdapter):
    def __init__(self, env="arc-stage2-v1-stub"):
        super().__init__("stub", env, ENTRY)
