"""Step 09 renders the shared three-metric Stage 2 publication."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from stage2.evaluation.figures import render
