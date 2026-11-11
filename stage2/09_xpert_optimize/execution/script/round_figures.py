"""Step 09 renders the shared three-metric Stage 2 publication."""
from pathlib import Path
import sys
sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())))
from stage2.evaluation.figures import render
