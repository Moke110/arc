"""Exercise the public fit and fresh-process replay contracts."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / 'output/api'))

from tsg_transfer.fit import replay_fit, train, verify_replay


def fixture(root):
    data = SimpleNamespace(directory=root, source='tahoe_bulk', fingerprints=lambda rows:
                           np.ones((len(rows), 2048), dtype=np.float32))
    rng = np.random.default_rng(42)
    data.expression = rng.uniform(size=(8, 5000)).astype(np.float32)
    data.control = rng.uniform(size=(8, 5000)).astype(np.float32)
    data.doses = lambda rows: np.asarray([.05, .5, 5., .5, .05, .5, 5., .5])[rows]
    return data


class FitTest(unittest.TestCase):
    def test_bounded_native_fit_replays_the_best_checkpoint_in_a_fresh_process(self):
        workspace = Path('/home/chang/arc')
        source = workspace / 'pert_model/TranSiGen/source'
        if not source.exists():
            self.skipTest('native TranSiGen source unavailable')
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            (root / 'manifest.json').write_text('{}')
            candidate = root / 'fold-0'
            train_rows, validation_rows = np.arange(4), np.arange(4, 8)
            data = fixture(root)
            receipt = train(data, source, candidate, train_rows=train_rows,
                            validation_rows=validation_rows, fold=0, kind='tahoe',
                            device='cpu', smoke_epochs=2)
            self.assertGreater(receipt['training_steps'], 0)
            self.assertEqual(receipt['learning_rate'], 1e-4)
            self.assertLessEqual(receipt['best_epoch'], 2)
            self.assertAlmostEqual(receipt['transform']['mean'], np.log10(.5))
            with self.assertRaisesRegex(ValueError, 'fresh process'):
                replay_fit(data, source, candidate, train_rows=train_rows,
                           validation_rows=validation_rows)
            result = subprocess.run([sys.executable, __file__, str(root), str(source)],
                                    check=True, capture_output=True, text=True)
            proof = json.loads(result.stdout)
            self.assertEqual(receipt['validation_mse'], proof['validation_mse'])
            self.assertEqual(receipt, verify_replay(data, candidate, train_rows=train_rows,
                             validation_rows=validation_rows, kind='tahoe', fold=0, formal=False))
            with self.assertRaisesRegex(ValueError, 'formal'):
                verify_replay(data, candidate, train_rows=train_rows,
                              validation_rows=validation_rows, kind='tahoe', fold=0)
            (candidate / 'learning_curve.csv').write_text('tampered')
            with self.assertRaisesRegex(ValueError, 'digest'):
                verify_replay(data, candidate, train_rows=train_rows,
                              validation_rows=validation_rows, kind='tahoe', fold=0, formal=False)


if __name__ == '__main__':
    root, source = map(Path, sys.argv[1:])
    print(json.dumps(replay_fit(fixture(root), source, root / 'fold-0',
                                train_rows=np.arange(4), validation_rows=np.arange(4, 8))))
