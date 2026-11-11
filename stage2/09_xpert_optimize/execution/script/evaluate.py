"""Step 09 uses the shared Stage 2 full-CM significance scorer."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from stage2.evaluation.evaluate import main

if __name__ == '__main__':
    sys.argv.extend(['--step', '9'])
    main()
