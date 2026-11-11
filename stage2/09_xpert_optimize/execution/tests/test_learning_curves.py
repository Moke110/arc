"""The queue refreshes learning curves in the requested output on resumed runs."""
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
import round1_fit


class LearningCurvesTest(unittest.TestCase):
    def test_resumed_queue_refreshes_custom_output(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            work, output = root / 'work', root / 'published'
            policy = {'phases': {phase: {'updates': 20} for phase in ('tahoe', 'tenv')}}
            entries = [{'name': f'XPert_{i}', 'predictor_sha256': 'input'} for i in range(9)]
            plan = {'schema': 'step09-round1-plan-v1', 'fold': 0, 'seed': 17,
                    'models': entries, 'policy_sha256': 'policy',
                    'gene_axes_manifest_sha256': 'axis', 'fit_sources': {}}
            for entry in entries:
                fit = work / 'fits' / entry['name']
                fit.mkdir(parents=True)
                history = [{'update': update, 'train_terms': [1, 1, 1, 1],
                            'validation_loss': loss} for update, loss in [(10, 1.0), (20, .8)]]
                (fit / 'run.json').write_text(json.dumps(dict(
                    status='complete', model=entry['name'], policy=policy, policy_sha256='policy',
                    predictor_sha256='input', fold=0, seed=17,
                    phases={phase: {'history': history} for phase in ('tahoe', 'tenv')})))
            plan_path, policy_path = root / 'plan.json', root / 'policy.json'
            plan_path.write_text(json.dumps(plan))
            policy_path.write_text(json.dumps(policy))
            args = ['queue', '--role', 'validation', '--work-root', str(work), '--output', str(output),
                    '--plan', str(plan_path), '--policy', str(policy_path)]
            with patch.object(sys, 'argv', args), patch.object(round1_fit, 'verify_committed'), \
                    patch.object(round1_fit, 'check_digest'), patch.object(round1_fit, 'run_queue', return_value=[]):
                round1_fit.main()
                figure = output / 'xpert_opti_r1_learning_curves.png'
                first = figure.read_bytes()
                figure.write_bytes(b'stale plot')
                round1_fit.main()
            self.assertEqual(figure.read_bytes(), first)
            with (output / 'learning_curve_points.csv').open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 36)
            self.assertEqual({row['update'] for row in rows}, {'10', '20'})
            self.assertTrue((output / 'learning-curves.json').is_file())
            self.assertFalse(list(output.glob('xpert_opti_r1_learning_XPert_*.png')))
