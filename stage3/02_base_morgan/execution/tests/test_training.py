"""Synthetic saved-run acceptance at the public training helpers."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow import digest, encoded


class SavedRunTest(unittest.TestCase):
    def test_saved_checkpoint_inference_and_corruption(self):
        from training import train, verify, evaluate, predict
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = [dict(measurement_id=str(i), canonical_smiles='C' * (i + 1),
                         parent_id=str(i), condition='baseline', phenotype='cell_viability',
                         dose_um=dose, direction=(i % 3)-1)
                    for i, dose in enumerate([1, 10, 1, 10, 100, 1000])]
            splits = [dict(setting='cold_drug', fold=0,
                           roles=['train']*4+['validation', 'test'])]
            import rdkit
            identity = dict(condition='baseline', phenotype='cell_viability', setting='cold_drug', fold=0, seed=17)
            prepared = dict(rows=rows, splits=splits, input_fingerprint=digest(rows),
                            split_fingerprint=digest(splits), expected_runs=[identity],
                            features=dict(radius=2, bits=2048, chirality=True, rdkit_version=rdkit.__version__))
            source = root/'prepared.json'
            source.write_bytes(encoded(prepared))
            result = train(source, root/'run', **identity, smoke_epochs=2)
            self.assertEqual(result['dose_statistics']['mean'], 0.5)
            self.assertEqual(result['dose_statistics']['std'], 0.5)
            self.assertEqual(result['mode'], 'smoke')
            self.assertTrue(verify(root/'run', source))
            self.assertEqual(evaluate(root/'run')['n'], 1)
            probabilities = predict(root/'run'/'checkpoint.pt', [rows[-1]])
            self.assertAlmostEqual(sum(probabilities[0]), 1, places=6)
            records = json.loads((root/'run'/'predictions.json').read_text())
            self.assertEqual(probabilities, [records[0]['probabilities']])
            records[0]['direction'] = 99
            (root/'run'/'predictions.json').write_bytes(encoded(records))
            with self.assertRaisesRegex(ValueError, 'hash|fingerprint'):
                verify(root/'run', source)

if __name__ == '__main__':
    unittest.main()
