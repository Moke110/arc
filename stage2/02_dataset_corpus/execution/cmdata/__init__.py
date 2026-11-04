"""Stage 2 CM corpus build package.

The package is deliberately stdlib + pandas/numpy/scanpy/rdkit/scipy only; it
never imports Stage 1 code.  `build_cm_datasets.py` is the single entry point.
"""

__all__ = ["paths", "text", "annotation", "genes", "sources", "drugs",
           "views", "qa"]
