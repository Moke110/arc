#!/usr/bin/env python
"""Run PRnet's upstream test_demo.py on CPU without editing the clone.

The public LINCS checkpoint was saved from CUDA tensors; upstream calls
torch.load() without map_location and therefore fails on a CPU-only torch.
This wrapper only adds that map_location and then executes the upstream script.
"""
import runpy
import sys
from pathlib import Path

import torch

_orig_load = torch.load


def _cpu_load(*args, **kwargs):
    kwargs.setdefault("map_location", "cpu")
    return _orig_load(*args, **kwargs)


torch.load = _cpu_load

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2] / "pert_model/prnet"
sys.path.insert(0, str(REPO))
sys.argv = [str(REPO / "test_demo.py"), *sys.argv[1:]]
runpy.run_path(str(REPO / "test_demo.py"), run_name="__main__")
