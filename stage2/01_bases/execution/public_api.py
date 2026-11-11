"""Bind the single public implementation published under ``output/api``.

ADR 0062 and the Stage 2 step contract place public reader/metric code in the
output contract; build, audit and test commands import that one implementation
instead of keeping a second copy beside them. Importing this module puts
``output/api`` on ``sys.path`` so ``import bases...`` resolves to the published
modules. Nothing here computes or caches scientific results.
"""
import sys
from pathlib import Path

API_ROOT = Path(__file__).resolve().parents[1] / "output" / "api"
if str(API_ROOT) not in sys.path:
    sys.path.insert(0, str(API_ROOT))
