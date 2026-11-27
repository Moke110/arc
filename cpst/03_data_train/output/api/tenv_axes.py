"""Resolve active tenv axes through the common plate-pooled training reader."""

from pathlib import Path

from training_data import DEFAULT, TENV_AXIS_STEPS as ORIGINAL_TENV_AXIS_STEPS, open_tenv_training

NEW_TENV_AXIS_STEPS = {
    "hvg-2k": "03.8_tenv_hvg-2k",
    "denoise-2k": "03.9_tenv_denoise-2k",
}
TENV_AXIS_STEPS = ORIGINAL_TENV_AXIS_STEPS | NEW_TENV_AXIS_STEPS


def tenv_axis_output(axis_name, root=None):
    """Resolve a named active tenv release under a Step 03 root."""
    if axis_name not in TENV_AXIS_STEPS:
        raise ValueError(f"unknown tenv training axis: {axis_name}")
    step = DEFAULT.parent if root is None else Path(root)
    return step / TENV_AXIS_STEPS[axis_name] / "output"


def open_tenv_axis_training(axis_name, root=None):
    """Verify and open one named active tenv axis under a Step 03 root."""
    return open_tenv_training(tenv_axis_output(axis_name, root))
