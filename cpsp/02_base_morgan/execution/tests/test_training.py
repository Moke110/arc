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
            from rdkit import Chem
            from workflow import cold_drug, SPLIT_SEED, TRAINING_SEEDS, CLASSES
            from training import _feature_configuration
            for row in rows:
                row['parent_id'] = Chem.MolToInchiKey(Chem.MolFromSmiles(row['canonical_smiles']))
            splits = cold_drug(rows)
            # Shift held-out doses without changing parent membership.
            for i, role in enumerate(splits[0]['roles']):
                rows[i]['dose_um'] = 1000 if role == 'test' else 100 if role == 'validation' else (1 if i % 2 else 10)
            splits = cold_drug(rows)
            identity = dict(condition='baseline', phenotype='cell_viability', setting='cold_drug', fold=0, seed=17)
            prepared = dict(schema_version=1, scope='cold_drug', split_seed=SPLIT_SEED,
                            training_seeds=TRAINING_SEEDS, class_order=CLASSES,
                            rows=rows, splits=splits, input_fingerprint=digest(rows),
                            split_fingerprint=digest(splits),
                            expected_runs=[dict(identity, fold=f, seed=seed) for f in range(3) for seed in TRAINING_SEEDS],
                            features=_feature_configuration())
            source = root/'prepared.json'
            source.write_bytes(encoded(prepared))
            import copy
            from training import _prepared
            corrupt = copy.deepcopy(prepared)
            corrupt['splits'][0]['roles'][0] = 'train' if splits[0]['roles'][0] != 'train' else 'test'
            corrupt['split_fingerprint'] = digest(corrupt['splits'])
            source.write_bytes(encoded(corrupt))
            with self.assertRaisesRegex(ValueError, 'canonical'):
                train(source, root/'rejected-run', **identity, smoke_epochs=1)
            self.assertFalse((root/'rejected-run').exists())
            corrupt = copy.deepcopy(prepared)
            del corrupt['features']['effective_configuration']
            source.write_bytes(encoded(corrupt))
            with self.assertRaisesRegex(ValueError, 'configuration'):
                _prepared(source)
            source.write_bytes(encoded(prepared))
            result = train(source, root/'run', **identity, smoke_epochs=2)
            from training import _statistics
            training_rows = [r for r, role in zip(rows, splits[0]['roles']) if role == 'train']
            self.assertEqual(result['dose_statistics'], _statistics(training_rows))
            self.assertEqual(result['mode'], 'smoke')
            self.assertTrue(verify(root/'run', source))
            self.assertEqual(evaluate(root/'run')['n'], 2)
            probabilities = predict(root/'run'/'checkpoint.pt', [r for r, role in zip(rows, splits[0]['roles']) if role == 'test'])
            self.assertAlmostEqual(sum(probabilities[0]), 1, places=6)
            records = json.loads((root/'run'/'predictions.json').read_text())
            self.assertEqual(probabilities, [r['probabilities'] for r in records])
            records[0]['direction'] = 99
            (root/'run'/'predictions.json').write_bytes(encoded(records))
            with self.assertRaisesRegex(ValueError, 'hash|fingerprint'):
                verify(root/'run', source)

if __name__ == '__main__':
    unittest.main()
