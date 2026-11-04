"""Stage 2 dataset-corpus build package.

The package is deliberately stdlib + pandas/numpy/scipy/rdkit/anndata only; it
never imports Stage 1 code.  `build_corpus.py` is the single entry point; the
four CM corpora are written by `corpus_h5ad.py` (native gene axes, ADR 0033).
"""

__all__ = ["paths", "text", "annotation", "genes", "sources", "drugs",
           "corpus", "corpus_h5ad", "qa"]
