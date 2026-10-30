#!/usr/bin/env python
"""Run PRnet's upstream test_demo.py on CPU without editing the clone.

The public LINCS checkpoint was saved from CUDA tensors; upstream calls
torch.load() without map_location and therefore fails on a CPU-only torch.
This wrapper only adds that map_location and then executes the upstream script.
"""

import sys

# Keep `pert_model/*/source/` byte-identical to its download (ADR 0030): this
# smoke imports model code from the vendored tree, so byte-code must not land
# there. Set before any other import; the run.sh wrapper exports the same flag.
sys.dont_write_bytecode = True

import runpy
from pathlib import Path

import torch

_orig_load = torch.load


def _cpu_load(*args, **kwargs):
    kwargs.setdefault("map_location", "cpu")
    return _orig_load(*args, **kwargs)


torch.load = _cpu_load

HERE = Path(__file__).resolve().parent


def find_repo_root(start: Path) -> Path:
    """Walk up to the directory containing .git (ADR 0030 §4)."""
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RuntimeError(f"cannot find the repository root above {start}")


REPO = find_repo_root(HERE) / "pert_model/PRnet/source"
sys.path.insert(0, str(REPO))
sys.argv = [str(REPO / "test_demo.py"), *sys.argv[1:]]
runpy.run_path(str(REPO / "test_demo.py"), run_name="__main__")
