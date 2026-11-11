"""Formal queue refuses partial, fixture and altered completed runs."""
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from rerun_queue import MODELS, verify_plan, completed
from round_contract import sha256


class QueueTest(unittest.TestCase):
    def test_partial_or_fixture_plan_cannot_launch_formal_training(self):
        plan = {'schema': 'step09-round1-rerun-plan-v1', 'fixture': False, 'runs': []}
        with self.assertRaises(ValueError):
            verify_plan(plan)
        plan['runs'] = [{'id': f'{m}-{f}-{s}', 'model': m, 'fold': f, 'seed': s}
                        for m in MODELS for f in (0, 1) for s in (17, 42)]
        plan['fixture'] = True
        with self.assertRaises(ValueError):
            verify_plan(plan)

    def test_completed_fit_requires_selected_checkpoint_identity(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            policy = root / 'policy.json'
            policy.write_text('{}')
            run = {'id': 'x', 'model': 'XPert_G512', 'fold': 0, 'seed': 17, 'predictor_sha256': 'p'}
            (root / 'best.pt').write_bytes(b'selected')
            report = {'status': 'complete', **run, 'policy_sha256': sha256(policy),
                      'phases': {'tenv': {'checkpoint_sha256': sha256(root / 'best.pt')}}}
            (root / 'run.json').write_text(json.dumps(report))
            self.assertEqual(completed(run, root, policy)['status'], 'complete')
            (root / 'best.pt').write_bytes(b'changed')
            with self.assertRaises(ValueError):
                completed(run, root, policy)
