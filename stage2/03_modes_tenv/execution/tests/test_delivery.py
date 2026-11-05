"""Publication through the artifact-only CLI; no model fitting or test truth."""
import json
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from artifacts import sha256, write_csv, write_json
from metrics import METRICS

SCRIPT = Path(__file__).resolve().parents[1] / 'script' / 'deliver.py'


def evaluation_fixture(output):
    output.mkdir()
    summaries = []
    runs = []
    bands = []
    selections = []
    for model, completed in [('no_change', 9), ('delta_mean', 8), ('failed_model', 0)]:
        summary = {'model': model, 'completed_runs': completed, 'expected_runs': 9,
                   'complete_ranking': completed == 9, 'deterministic_reference': model != 'failed_model'}
        for metric in METRICS:
            summary.update({f'{metric}_mean': .1 if completed else np.nan,
                            f'{metric}_sample_sd': .01 if completed else np.nan,
                            f'{metric}_valid_runs': completed})
        summaries.append(summary)
        for index in range(9):
            runs.append({'model': model, 'fold': index // 3, 'seed': (17, 42, 2026)[index % 3],
                         'status': 'completed' if index < completed else 'timeout',
                         **{metric: .1 if index < completed else np.nan for metric in METRICS}})
        for direction in ('positive', 'negative'):
            for recall in (0., .5, 1.):
                bands.append({'model': model, 'direction': direction, 'recall': recall,
                              'mean': .1 if completed else np.nan, 'min': .05 if completed else np.nan,
                              'max': .15 if completed else np.nan, 'valid_runs': completed})
        selections.append({'model': model, 'fold': 0 if completed else '',
                           'seed': 17 if completed else '', 'selection_score': .1 if completed else np.nan,
                           'label': 'post-hoc best test; fixed fold/seed tie order'})
    write_csv(output / 'summary.csv', summaries)
    write_csv(output / 'per_run.csv', runs)
    write_csv(output / 'per_group.csv', [{'model': 'no_change', 'group_id': 'g', 'PCC_reason': 'constant prediction'}])
    write_csv(output / 'best_runs.csv', selections)
    write_csv(output / 'pr_bands.csv', bands)
    np.savez(output / 'pr_sources.npz', recall=np.array([0., .5, 1.]),
             identifiers=np.array(['no_change:0:17', 'delta_mean:0:17']),
             curves=np.array([[[.1, .1, .1], [.1, .1, .1]], [[.2, .1, .1], [.1, .1, .1]]]))
    write_json(output / 'predictions.json', {'sealed': True})
    names = ['summary.csv', 'per_run.csv', 'per_group.csv', 'best_runs.csv', 'pr_bands.csv',
             'pr_sources.npz', 'predictions.json']
    write_json(output / 'evaluation.json', {'artifacts': {name: sha256(output / name) for name in names}})


class Delivery(unittest.TestCase):
    def test_publish_preserves_scored_artifacts_and_labels_best_run_limits(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'output'
            publication = root / 'publication'
            evaluation_fixture(output)
            original = {p.name: sha256(p) for p in output.iterdir()}
            result = subprocess.run([sys.executable, str(SCRIPT), '--output', str(output),
                                     '--publish', str(publication)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(original, {name: sha256(output / name) for name in original})
            text = ' '.join(''.join(node.itertext()) for node in ET.parse(publication / 'best_run.svg').iter()
                            if node.tag.endswith('}text'))
            self.assertIn('no_change (completed 9/9) [deterministic copies]', text)
            self.assertIn('delta_mean (completed 8/9; incomplete) [deterministic copies]', text)
            self.assertIn('failed_model (completed 0/9; incomplete): undefined selection', text)
            self.assertIn('Post-hoc best test run', text)
            manifest = json.loads((publication / 'publication.json').read_text())
            self.assertEqual(len(manifest['artifacts']), 11)
            for name, digest in manifest['artifacts'].items():
                self.assertEqual(sha256(publication / name), digest)
                self.assertEqual(sha256(output / name), digest)
            self.assertEqual(manifest['evaluation_sha256'], original['evaluation.json'])
            self.assertEqual(manifest['source_artifacts'], json.loads((output / 'evaluation.json').read_text())['artifacts'])

    def test_changed_source_is_rejected_before_publication(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'output'
            publication = root / 'publication'
            evaluation_fixture(output)
            (output / 'summary.csv').write_text('tampered table\n')
            result = subprocess.run([sys.executable, str(SCRIPT), '--output', str(output),
                                     '--publish', str(publication)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('replay checksum mismatch: summary.csv', result.stderr)
            self.assertFalse(publication.exists())
            self.assertFalse((output / 'metrics.png').exists())


if __name__ == '__main__':
    unittest.main()
