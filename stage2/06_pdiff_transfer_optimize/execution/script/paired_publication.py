"""Step 06 reporting uses the shared Stage 2 significance evaluator."""
from pathlib import Path
import sys
sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())))
from stage2.evaluation.evaluate import main

if __name__ == '__main__':
    sys.argv.extend(['--step', '6'])
    main()
