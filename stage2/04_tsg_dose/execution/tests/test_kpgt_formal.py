"""Public tenv KPGT formal planner contract."""

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

DRIVER = Path(__file__).resolve().parents[1] / 'kpgt/formal.py'


class FormalPlannerTest(unittest.TestCase):
    def test_plans_all_tenv_runs_and_keeps_the_baselines_separate(self):
        run = subprocess.run(
            [sys.executable, str(DRIVER), 'plan'],
            capture_output=True,
            text=True,
            check=True,
        )
        plans = json.loads(run.stdout)
        methods = {
            'no-dose',
            'scaling',
            'input-numeric',
            'latent-numeric',
            'input-embedding',
            'latent-embedding',
            'film',
        }
        expected = {
            (scenario, fold, method, 'tenv', 42)
            for scenario in ('cold-drug', 'cold-drug-dose')
            for fold in range(3)
            for method in methods
        }

        def coordinates(items):
            return {
                (
                    item['scenario'],
                    item['fold'],
                    item['method'],
                    item['source'],
                    item['seed'],
                )
                for item in items
            }

        self.assertEqual(len(plans['transigen']), 42)
        self.assertEqual(coordinates(plans['transigen']), expected)
        self.assertEqual(len(plans['baselines']), 12)
        self.assertEqual(
            coordinates(plans['baselines']),
            {
                (scenario, fold, method, 'tenv', 42)
                for scenario in ('cold-drug', 'cold-drug-dose')
                for fold in range(3)
                for method in ('KPGT_MLP', 'Mean_cell')
            },
        )


class FigureReplayTest(unittest.TestCase):
    def test_rejects_incomplete_run_grid(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            for table in ('pooled', 'folds', 'learning', 'coverage'):
                (root / f'{table}.csv').write_text('scenario,fold,method\n')
            runner = DRIVER.with_name('render.py')
            result = subprocess.run(
                [
                    sys.executable,
                    str(runner),
                    '--data',
                    str(root),
                    '--output',
                    str(root / 'figures'),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                'pooled.csv has an incomplete or duplicate plot grid', result.stderr
            )
            self.assertFalse((root / 'figures').exists())


class PublicationAcceptanceTest(unittest.TestCase):
    def test_accept_rejects_formal_plan_changed_after_preparation(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = root / 'inputs'
            formal = root / 'formal'
            inputs.mkdir()
            formal.mkdir()
            manifest = inputs / 'manifest.json'
            plan = formal / 'formal_plan.json'
            manifest.write_bytes(b'frozen input identity')
            plan.write_bytes(b'locked formal plan')
            receipt = root / 'candidate.json'
            receipt.write_text(
                json.dumps(
                    dict(
                        schema='tsg-dose-tenv-kpgt-report-v1',
                        status='candidate-verified',
                        input_manifest_sha256=sha256(manifest.read_bytes()).hexdigest(),
                        formal_plan_sha256=sha256(plan.read_bytes()).hexdigest(),
                    )
                )
            )
            plan.write_bytes(b'changed formal plan')
            acceptance = root / 'acceptance.json'
            result = subprocess.run(
                [
                    sys.executable,
                    str(DRIVER.with_name('publish.py')),
                    'accept',
                    '--candidate',
                    str(root / 'report'),
                    '--receipt',
                    str(receipt),
                    '--acceptance',
                    str(acceptance),
                    '--visual-inspection',
                    '--inputs',
                    str(inputs),
                    '--formal',
                    str(formal),
                    '--source',
                    str(root / 'source'),
                    '--baselines',
                    str(root / 'baselines'),
                    '--resources',
                    str(root / 'resources'),
                    '--monitor',
                    str(root / 'monitor.log'),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                'sealed input generation or formal plan changed after preparation',
                result.stderr,
            )
            self.assertFalse(acceptance.exists())


if __name__ == '__main__':
    unittest.main()
