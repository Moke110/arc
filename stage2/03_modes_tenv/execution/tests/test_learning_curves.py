"""Lightweight history mapping contracts; no model/data/evaluator imports."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('learning_curves', Path(__file__).parents[1] / 'script/learning_curves.py')
curves = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(curves)


class HistoryMappingTests(unittest.TestCase):
    def test_native_keys_and_actual_coordinates(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            for model, (x, train, validation, _, _) in curves.MAPPINGS.items():
                directory = root / 'runs' / model / 'fold_0/seed_17'
                directory.mkdir(parents=True)
                (directory / 'checkpoint.json').write_text(json.dumps({
                    'identity': {'model': model, 'fold': 0, 'seed': 17},
                    'history': [{x: 98, train: -2, validation: -3},
                                {x: 203, train: -4, validation: -5}]}))
            records = curves.collect(root)
            self.assertEqual(len(records), 63)
            available = [r for r in records if r['history_status'] == 'available']
            self.assertEqual(len(available), 5)
            for record in available:
                self.assertEqual(record['train'], [(98., -2.), (203., -4.)])
                self.assertEqual(record['validation'], [(98., -3.), (203., -5.)])
            self.assertEqual(sum(r['history_status'] == 'missing' for r in records), 40)
            self.assertEqual(sum(r['history_status'] == 'not applicable' for r in records), 18)

    def test_no_fabricated_coordinates_or_nonfinite_values(self):
        for history in ([{'step': 3, 'loss': 1}, {'step': 2, 'loss': 1}],
                        [{'step': 3, 'loss': float('nan')}]):
            with self.assertRaises(ValueError):
                curves.points(history, 'step', 'loss')
        self.assertEqual(curves.points([], 'epoch', 'loss'), [])

    def test_checkpoint_identity_is_required(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            directory = root / 'runs/prnet/fold_0/seed_17'
            directory.mkdir(parents=True)
            (directory / 'checkpoint.json').write_text('{"identity": {"model": "cpa"}}')
            with self.assertRaisesRegex(ValueError, 'identity mismatch'):
                curves.collect(root)


if __name__ == '__main__':
    unittest.main()
