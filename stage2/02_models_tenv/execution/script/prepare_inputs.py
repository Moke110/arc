"""CLI entry point for the bounded predictor-input build (#163).

The shared environment (``/home/chang/miniconda3/envs/arc/bin/python``) runs
this wrapper because the accepted ``bases`` readers need its newer
``anndata``.  The pinned PRnet environment never reads the corpus.
"""

import sys

from inputs import main

if __name__ == '__main__':
    sys.exit(main())
