"""Partial artifact boundary: real scoring, immutable source, no final seal."""
import sys
import tempfile
import unittest
from pathlib import Path
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from test_pipeline import fixture
from artifacts import SEEDS, read_json, write_json, rows, sha256, run_identity, run_directory
from benchmark import run_one, verify_barrier
from partial_benchmark import evaluate_partial, snapshot, verify_snapshot
from partial_deliver import deliver_partial


class PartialEvaluation(unittest.TestCase):
    def prepare(self, root):
        inputs, source = root / 'input', root / 'source'
        fixture(inputs)
        models = ['no_change', 'delta_mean', 'missing1', 'missing2', 'missing3', 'missing4', 'missing5']
        expected = [{'model': m, 'fold': f, 'seed': s,
                     'identity': run_identity(inputs, m, f, s)}
                    for m in models for f in range(3) for s in SEEDS]
        write_json(root / 'expected.json', expected)
        for model in models[:2]:
            run_one(inputs, source, model, 0, 17)
        return inputs, source, root / 'expected.json'

    def test_partial_values_counts_figures_source_preservation_and_seal_rejection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs, source, expected = self.prepare(root)
            before = {str(p): sha256(p) for p in source.rglob('*') if p.is_file()}
            output = root / 'partial'
            summary = evaluate_partial(inputs, source, output, expected, {}, 1)
            self.assertEqual(summary[0]['completed_runs'], 1)
            self.assertFalse(summary[0]['complete_ranking'])
            self.assertEqual(summary[2]['completed_runs'], 0)
            self.assertEqual(float(rows(output / 'per_run.csv')[0]['AUPRC_positive@100']), .5)
            self.assertEqual(int(rows(output / 'per_run.csv')[0]['PCC_valid_groups']), 0)
            self.assertEqual(len(rows(output / 'per_group.csv')), 2)
            self.assertEqual(before, {str(p): sha256(p) for p in source.rglob('*') if p.is_file()})
            self.assertFalse(read_json(output / 'predictions.json')['sealed'])
            with self.assertRaisesRegex(ValueError, 'barrier'):
                verify_barrier(inputs, output)
            for name in ('metrics', 'pr_panels', 'best_run'):
                texts = ' '.join(''.join(n.itertext()) for n in ET.parse(output / (name + '.svg')).iter()
                                 if n.tag.endswith('}text'))
                self.assertIn('PARTIAL finished-run snapshot', texts)
                self.assertGreater((output / (name + '.png')).stat().st_size, 1000)
            publication = deliver_partial(output, source.parent)
            self.assertEqual(len(publication['artifacts']), 6)
            for name, digest in publication['artifacts'].items():
                self.assertEqual(sha256(output / name), digest)
                self.assertEqual(sha256(root / name), digest)
            texts = ' '.join(''.join(n.itertext()) for n in ET.parse(root / 'best_run.svg').iter()
                             if n.tag.endswith('}text'))
            self.assertIn('available completed runs only', texts)
            self.assertIn('1/9 completed; incomplete, descriptive only', texts)
            self.assertIn('Nine-run entries: primary', texts)
            self.assertIn('completed runs:', texts)
            with self.assertRaisesRegex(ValueError, 'fixed step output'):
                deliver_partial(output, root / 'other')
            with (output / 'per_run.csv').open('a') as stream:
                stream.write('tampered')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                deliver_partial(output, root)
            with self.assertRaisesRegex(ValueError, 'fresh'):
                evaluate_partial(inputs, source, output, expected, {}, 1)
            with self.assertRaisesRegex(ValueError, 'separate'):
                evaluate_partial(inputs, source, source / 'partial', expected, {}, 1)

    def test_frozen_identity_and_changed_completed_prediction_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs, source, expected = self.prepare(root)
            captured = snapshot(inputs, source, expected, {})
            prediction = run_directory(source, 'no_change', 0, 17) / 'prediction.npy'
            with prediction.open('ab') as stream:
                stream.write(b'tampered')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                verify_snapshot(inputs, source, captured)
            frozen = read_json(expected)
            frozen[0]['identity']['input_sha256'] = 'wrong'
            write_json(expected, frozen)
            with self.assertRaisesRegex(ValueError, 'frozen expected'):
                snapshot(inputs, source, expected, {})


if __name__ == '__main__':
    unittest.main()
