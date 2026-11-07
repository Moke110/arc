"""Tahoe's frozen pilot protocol through the public planning CLI."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

EXECUTION = Path(__file__).resolve().parents[1]
CLI = EXECUTION / 'script/pilot.py'
PROTOCOL = EXECUTION / 'tahoe_pilot_record/protocol.json'


class TahoePilotTest(unittest.TestCase):
    def test_plan_uses_frozen_tahoe_protocol_and_schedules_all_fold_requests(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'pilot'
            release = root / 'verified-release'
            trainer = EXECUTION / 'script/train.py'
            result = subprocess.run(
                [sys.executable, str(CLI), 'profile', str(output), '--tahoe-release', str(release),
                 '--protocol', str(PROTOCOL), '--trainer', str(trainer)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            plan = json.loads((output / 'plan.json').read_text())
            self.assertEqual(plan['dataset'], 'tahoe')
            self.assertEqual(plan['protocol'], json.loads(PROTOCOL.read_text()))
            queue = json.loads((output / 'profile/queue.json').read_text())
            preparations = queue['jobs'][:3]
            self.assertEqual([job['kind'] for job in preparations], ['cpu'] * 3)
            for fold, prep in enumerate(preparations):
                request = plan['requests'][str(fold)]
                self.assertEqual(prep['command'][1], str(trainer))
                self.assertEqual(prep['command'][-5:],
                                 ['prepare-tahoe', str(release), request, '--fold', str(fold)])
                jobs = [job for job in queue['jobs'][3:] if job['id'].startswith(f'fold-{fold}-')]
                self.assertEqual(len(jobs), 2)
                recipes = []
                for job in jobs:
                    self.assertEqual(job['command'][1], str(trainer))
                    self.assertEqual(job['command'][-4:-2], ['fit', request])
                    recipes.append(json.loads(Path(job['command'][-2]).read_text()))
                self.assertEqual([recipe.pop('validation_interval') for recipe in recipes], [250, 500])
                self.assertEqual(recipes[0], recipes[1])
                self.assertEqual(recipes[0]['max_steps'], 1000)
                self.assertGreater(recipes[0]['minimum_steps'], 1000)


if __name__ == '__main__':
    unittest.main()
