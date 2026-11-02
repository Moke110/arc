"""Deterministic seeding and device resolution shared by the CLIs.

Seeding is part of the recorded run contract: python/numpy/torch RNGs are seeded
with the run seed, cuDNN is put in deterministic mode, and the training noise is
drawn from a CPU generator that is itself seeded from the run seed (see
``ptlib.train``).  ``resolve_device`` refuses a silent CPU fallback when CUDA was
requested: a CPU run is a different claim and must be asked for explicitly.
"""

from __future__ import annotations

import logging
import platform
import random

import numpy as np
import torch

__all__ = ["seed_everything", "resolve_device", "setup_logging"]


def seed_everything(seed: int) -> None:
    random.seed(int(seed))
    np.random.seed(int(seed) % (2**32 - 1))
    torch.manual_seed(int(seed))
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(int(seed))
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def resolve_device(requested: str) -> tuple[torch.device, str]:
    available = torch.cuda.is_available()
    if str(requested).startswith("cuda"):
        if not available:
            raise SystemExit(
                "--device cuda requested but CUDA is unavailable: refusing a silent CPU "
                "fallback (a CPU run is not the pinned GPU evidence). Use --device cpu "
                "explicitly for a CPU run."
            )
        index = 0
        if ":" in str(requested):
            index = int(str(requested).split(":", 1)[1])
        return torch.device(f"cuda:{index}"), torch.cuda.get_device_name(index)
    return torch.device("cpu"), platform.processor() or "cpu"


def setup_logging(run_log: Path | str | None, level: int = logging.INFO) -> logging.Logger:
    """Log to stderr and (when given) one run log file; overwrite, never append."""
    import logging
    import sys
    from pathlib import Path

    root = logging.getLogger()
    root.setLevel(level)
    for handler in list(root.handlers):
        root.removeHandler(handler)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    stream = logging.StreamHandler(sys.stderr)
    stream.setFormatter(formatter)
    root.addHandler(stream)
    if run_log is not None:
        path = Path(run_log)
        path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(path, mode="w", encoding="utf-8")
        file_handler.setFormatter(formatter)
        root.addHandler(file_handler)
    return root
