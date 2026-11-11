"""Stage 2 plots use only the shared three-metric significance renderer."""
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())))
from stage2.evaluation.figures import render
from stage2.evaluation.metrics import METRICS
