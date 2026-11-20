"""Stage 2 TranSiGen transfer, using the accepted Step 03 model and inputs."""
from pathlib import Path
import sys

STEP_03 = Path(__file__).resolve().parents[4] / '03_tsg_dose'
step_03_api = str(STEP_03 / 'output/api')
if step_03_api not in sys.path:
    sys.path.insert(0, step_03_api)
