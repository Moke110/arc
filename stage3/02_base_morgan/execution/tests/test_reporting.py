"""Known saved predictions exercise reporting without a training dependency."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reporting import _digest, summarize, render_report


class ReportingTests(unittest.TestCase):
    def fixture(self, root):
        rows = [dict(measurement_id=str(i), parent_id=str(i), dose_um=1, direction=-1,
                     condition='baseline', phenotype='cell_viability', source_cell=str(i)) for i in range(4)]
        splits = [dict(setting='cold_drug', fold=f, roles=['test' if (i if i < 2 else 2) == f else 'train' for i in range(4)]) for f in range(3)]
        identities = [dict(condition='baseline', phenotype='cell_viability', setting='cold_drug', fold=f, seed=17) for f in range(3)]
        prepared = dict(rows=rows, splits=splits, expected_runs=identities, input_fingerprint=_digest(rows), split_fingerprint=_digest(splits), features={})
        paths = []
        for identity, split in zip(identities, splits):
            path = root / str(identity['fold'])
            path.mkdir()
            records = [dict(r, run_identity=identity, probabilities=([1, 0, 0] if identity['fold'] < 2 else [0, 1, 0]), prediction=(-1 if identity['fold'] < 2 else 0)) for r, role in zip(rows, split['roles']) if role == 'test']
            (path/'predictions.json').write_text(json.dumps(records))
            manifest = dict(identity=identity, status='completed', mode='formal', smoke_epochs=None, features={}, input_fingerprint=prepared['input_fingerprint'], split_fingerprint=prepared['split_fingerprint'], artifact_hashes={'predictions.json': hashlib.sha256((path/'predictions.json').read_bytes()).hexdigest()})
            (path/'run.json').write_text(json.dumps(manifest))
            paths.append(path)
        return prepared, paths

    def test_absent_class_pooling_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            prepared, paths = self.fixture(Path(tmp))
            result = summarize(prepared, paths, require_complete=True)
            self.assertTrue(result['complete'])
            self.assertAlmostEqual(result['aggregates'][0]['macro_f1']['mean'], 2/9)
            self.assertAlmostEqual(result['cold_drug_pooled_oof'][0]['metrics']['macro_f1'], 2/9)
            # Loss weighs the two-error fold twice when measurements are pooled.
            self.assertNotAlmostEqual(result['aggregates'][0]['log_loss']['mean'], result['cold_drug_pooled_oof'][0]['metrics']['log_loss'])
            self.assertEqual(result['runs'][0]['metrics']['recall'], [1, 0, 0])
            json.dumps(result, allow_nan=False)
            self.assertIn('ddof=1', render_report(result))

    def test_missing_duplicate_unexpected_and_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            prepared, paths = self.fixture(Path(tmp))
            self.assertEqual(summarize(prepared, paths[:1])['status_counts']['pending'], 2)
            with self.assertRaisesRegex(ValueError, 'Incomplete'):
                summarize(prepared, paths[:1], require_complete=True)
            with self.assertRaisesRegex(ValueError, 'Duplicate run'):
                summarize(prepared, paths + paths[:1])
            manifest_path = paths[0]/'run.json'
            manifest = json.loads(manifest_path.read_text())
            manifest.update(mode='smoke', smoke_epochs=2)
            manifest_path.write_text(json.dumps(manifest))
            self.assertFalse(summarize(prepared, paths)['complete'])
            manifest['identity']['seed'] = 99
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'Unexpected'):
                summarize(prepared, paths)

    def test_readable_tables_and_mode_separation(self):
        with tempfile.TemporaryDirectory() as tmp:
            prepared, paths = self.fixture(Path(tmp))
            prepared['provenance_metadata'] = {'coverage': {'baseline': {
                'source_values': 7, 'admitted_values': 4, 'pending_values': 3, 'resolved_labels': 4}}}
            manifest_path = paths[0] / 'run.json'
            manifest = json.loads(manifest_path.read_text())
            manifest.update(mode='smoke', smoke_epochs=2)
            manifest_path.write_text(json.dumps(manifest))
            summary = summarize(prepared, paths[:2])
            held = dict(summary['aggregates'][0], setting='cold_dose', scope='held_dose',
                        held_dose=0.2, generalization_type='interpolation')
            summary['aggregates'].append(held)
            report = render_report(summary)
            self.assertIn('INCOMPLETE formal benchmark', report)
            self.assertIn('| 3 | 1 | 1 | 1 | 0 |', report)
            self.assertIn('| baseline | 7 | 4 | 3 | 4 |', report)
            self.assertIn('| baseline/cell_viability | cold_drug | formal |', report)
            self.assertIn('| baseline/cell_viability | cold_drug | smoke |', report)
            self.assertIn('| baseline/cell_viability | formal | 0.2 | interpolation |', report)
            self.assertLess(report.index('## Cold-dose held levels'), report.index('## Unweighted fold'))
            self.assertIn('## Per-class precision and recall', report)
            self.assertIn('## Cold-drug pooled OOF per seed', report)
            self.assertIn('[summary.json](summary.json)', report)
            self.assertIn('[run ledger](benchmark.json)', report)
            self.assertNotIn('```json', report)
            self.assertIn('not biological replication', report)
            self.assertIn('training-only', report)

    def test_corruption_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            prepared, paths = self.fixture(Path(tmp))
            path = paths[0]/'predictions.json'
            records = json.loads(path.read_text())
            records[0]['direction'] = 1
            path.write_text(json.dumps(records))
            with self.assertRaisesRegex(ValueError, 'hash'):
                summarize(prepared, paths)
            manifest_path = paths[0]/'run.json'
            manifest = json.loads(manifest_path.read_text())
            manifest['artifact_hashes']['predictions.json'] = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'provenance'):
                summarize(prepared, paths)


if __name__ == '__main__':
    unittest.main()
