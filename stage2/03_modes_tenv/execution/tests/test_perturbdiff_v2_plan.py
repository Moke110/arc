"""Offline migration planning through the approved artifact CLI seam."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'script/perturbdiff_v2_plan.py'


class MigrationPlanning(unittest.TestCase):
    def test_preserves_other_jobs_and_debits_historical_training(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            states = root / 'states'
            states.mkdir()
            jobs = [{'model': model, 'fold': fold, 'seed': seed, 'command': [model]}
                    for model in ('perturbdiff', 'cpa', 'prnet', 'chemcpa', 'cmonge', 'no_change', 'delta_mean')
                    for fold in range(3) for seed in (17, 42, 2026)]
            for job in jobs[:9]:
                key = f"perturbdiff_fold{job['fold']}_seed{job['seed']}"
                (states / (key + '.json')).write_text(json.dumps({'training_seconds': 12476,
                    'status': 'completed', 'identity': {'seed': job['seed']}}))
            queue = root / 'queue.json'
            queue.write_text(json.dumps({'jobs': jobs}))
            before = queue.read_bytes()
            output = root / 'plan.json'
            result = subprocess.run([sys.executable, str(SCRIPT), '--queue', str(queue),
                '--states', str(states), '--output', str(output)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            plan = json.loads(output.read_text())
            self.assertEqual(queue.read_bytes(), before)
            self.assertEqual(plan['untouched_jobs'], jobs[9:])
            self.assertEqual(len(plan['migrations']), 9)
            self.assertEqual(plan['migrations'][0]['remaining_training_seconds'], 1924)
            self.assertEqual(plan['migrations'][0]['initialization'], 'random_init')
            self.assertFalse(plan['dispatch_authorized'])
