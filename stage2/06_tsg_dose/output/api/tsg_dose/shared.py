"""Bind the accepted shared learner mechanisms without editing its namespace."""
from pathlib import Path
import sys
SHARED=Path(__file__).resolve().parents[4]/'04_models_tenv-5k'/'output/api'
if str(SHARED) not in sys.path:sys.path.insert(0,str(SHARED))
