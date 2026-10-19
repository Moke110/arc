"""
State adapter (ticket 11): drives the ST-HVG-Tahoe entry in its pinned env
`arc-stage2-v1-state` through the runner harness. All model access is through the
`predict` seam; this module never imports model internals.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from model_runner import SubprocessAdapter  # noqa: E402

ENTRY = Path(__file__).resolve().parent / "entry.py"


class StateAdapter(SubprocessAdapter):
    def __init__(self):
        super().__init__("state", "arc-stage2-v1-state", ENTRY)
