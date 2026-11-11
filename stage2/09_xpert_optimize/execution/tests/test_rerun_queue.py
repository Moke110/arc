"""Formal queue refuses partial, fixture and altered completed runs."""
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from rerun_queue import MODELS, verify_plan, completed, execute
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

    def test_only_two_models_form_the_complete_eight_run_plan(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'manifest.json').write_text('{}')
            runs = [{'id': f'{model}-{fold}-{seed}', 'model': model, 'fold': fold, 'seed': seed,
                     'predictor': str(root), 'predictor_sha256': sha256(root / 'manifest.json'),
                     'evaluator': str(root), 'evaluator_sha256': sha256(root / 'manifest.json')}
                    for model in ('XPert_G512', 'XPert_MLP512')
                    for fold in (0, 1) for seed in (17, 42)]
            plan = {'schema': 'step09-round1-rerun-plan-v1', 'fixture': False, 'runs': runs}
            verify_plan(plan)
            unwanted = dict(runs[0], model='XPert_G1024', id='excluded')
            with self.assertRaises(ValueError):
                verify_plan(dict(plan, runs=runs + [unwanted]))

    def test_persistent_user_cancellation_prevents_any_queue_launch(self):
        with TemporaryDirectory() as temp:
            work = Path(temp)
            (work / 'cancel-pending.json').write_text('{"status":"finish_current_then_stop"}')
            with self.assertRaises(SystemExit) as stopped:
                execute(work / 'unused-plan.json', work / 'unused-policy.json', work)
            self.assertEqual(stopped.exception.code, 75)
            self.assertEqual(json.loads((work / 'queue-status.json').read_text())['status'],
                             'cancelled_by_user')
            self.assertFalse((work / 'fits').exists())

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
