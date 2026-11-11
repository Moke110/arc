"""Retrieval reports use the shared scorer; exposure plotting is independent."""
from pathlib import Path
import sys
REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())
sys.path.insert(0, str(REPO))
from stage2.evaluation.evaluate import main
from exposure import render_exposure

if __name__ == '__main__':
    sys.argv.extend(['--step', '6'])
    main()
