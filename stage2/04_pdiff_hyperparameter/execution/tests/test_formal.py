"""Verify the formal producer rejects an altered plan before accessing inputs."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


CLI = Path(__file__).resolve().parents[1] / 'script/formal.py'
COLLECT = CLI.with_name('formal_outcomes.py')


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + '\n')
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FrozenWorker(unittest.TestCase):
    def test_worker_rejects_unpinned_plan_before_loading_native_runtime(self):
        with tempfile.TemporaryDirectory() as temporary:
            plan = Path(temporary) / 'plan.json'
            plan.write_text('{"altered": true}\n')
            for phase in ('fit', 'predict'):
                with self.subTest(phase=phase):
                    result = subprocess.run(
                        [sys.executable, str(CLI), phase, str(plan), '0' * 64, '0', '17'],
                        capture_output=True, text=True, timeout=30)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn('frozen plan hash mismatch', result.stderr)
                    self.assertEqual(list(plan.parent.iterdir()), [plan])


class TerminalOutcomes(unittest.TestCase):
    def test_reports_every_failed_fit_without_waiting_for_prediction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan_path = root / 'plan.json'
            folds = []
            for fold in range(3):
                request = f'fold-{fold}.json'
                folds.append({'request': request, 'request_sha256': write(root / request, {})})
            plan_sha = write(plan_path, {'schema': 'step04-evaluation-plan-v1', 'stage': 'formal',
                                         'dataset': 'tenv', 'folds': folds})
            for phase in ('fit', 'predict'):
                jobs = []
                for fold in range(3):
                    for seed in (17, 42, 2026):
                        jobs.append({'id': f'fold-{fold}-seed-{seed}', 'cwd': str(root),
                                     'kind': 'training' if phase == 'fit' else 'gpu',
                                     'memory_gib': 10, 'disk_gib': 2, 'timeout_seconds': 30,
                                     'command': [str(CLI), phase, str(plan_path), plan_sha, str(fold), str(seed)]})
                write(root / f'{phase}-queue.json', {'record_dir': f'{phase}-scheduler', 'jobs': jobs})
                if phase == 'fit':
                    record = root / 'fit-scheduler'
                    write(record / 'status.json', {'jobs': {
                        job['id']: {'identity': job, 'status': 'timeout', 'attempts': [
                            {'reason': 'worker timeout', 'log': str(record / (job['id'] + '.log'))}]}
                        for job in jobs}})
            output = root / 'outcomes.json'
            result = subprocess.run([sys.executable, str(COLLECT), str(plan_path), plan_sha, str(output)],
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            outcomes = json.loads(output.read_text())['outcomes']
            self.assertEqual(len(outcomes), 9)
            self.assertEqual({(row['phase'], row['status'], row['reason']) for row in outcomes},
                             {('fit', 'timeout', 'worker timeout')})
            self.assertTrue(all('prediction' not in row for row in outcomes))


class CheckpointRelease(unittest.TestCase):
    def test_publishes_all_nine_checkpoint_identities_and_replay_records(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            folds = []
            for fold in range(3):
                request_name = f'fold-{fold}/request.json'
                entry = {'request': request_name,
                         'request_sha256': write(root / request_name, {'identity': {'dataset': 'tenv', 'fold': fold}}),
                         'policy_sha256': 'b' * 64, 'recipe_sha256_by_seed': {}, 'recipes': {}}
                for seed in (17, 42, 2026):
                    name = f'fold-{fold}/seed-{seed}/recipe.json'
                    recipe = {'stage': 'formal', 'seed': seed, 'fold': fold}
                    entry['recipes'][str(seed)] = {'path': name, 'sha256': write(root / name, recipe)}
                    entry['recipe_sha256_by_seed'][str(seed)] = 'a' * 64
                folds.append(entry)
            plan = root / 'plan.json'
            plan_sha = write(plan, {'schema': 'step04-evaluation-plan-v1', 'stage': 'formal',
                                    'dataset': 'tenv', 'selection_sha256': 'e' * 64, 'folds': folds})
            rows = []
            for fold in range(3):
                for seed in (17, 42, 2026):
                    slot = root / f'fold-{fold}/seed-{seed}'
                    checkpoint = {'identity': {'dataset': 'tenv', 'fold': fold, 'seed': seed, 'stage': 'formal'},
                                  'checkpoint_sha256': write(slot / 'best.pt', {'fold': fold, 'seed': seed}),
                                  'recipe_sha256': 'a' * 64, 'best_step': 500,
                                  'best_validation_loss': 1., 'model_state_sha256': 'c' * 64}
                    reload = {'validation_loss': 1., 'prediction_max_abs_error': 0}
                    run_hash = write(slot / 'run.json', {**checkpoint, 'status': 'completed', 'reload': reload,
                                                        'upstream_commit': 'pinned', 'environment': {}})
                    replay_hash = write(slot / 'replay.json', {**checkpoint, **reload, 'plan_sha256': plan_sha})
                    prediction = slot / 'prediction.json'
                    prediction_hash = write(prediction, {**checkpoint, 'schema': 'step04-heldout-prediction-v1',
                        'plan_sha256': plan_sha, 'run_sha256': run_hash, 'replay_sha256': replay_hash})
                    rows.append({'dataset': 'tenv', 'fold': fold, 'seed': seed, 'status': 'completed',
                                 'prediction': str(prediction.relative_to(root)),
                                 'prediction_manifest_sha256': prediction_hash})
            outcomes = root / 'outcomes.json'
            outcomes_sha = write(outcomes, {'schema': 'step04-evaluation-outcomes-v1',
                                             'plan_sha256': plan_sha, 'outcomes': rows})
            evaluation = root / 'evaluation'
            write(evaluation / 'evaluation.json', {'schema': 'step04-evaluation-v1', 'stage': 'formal',
                  'complete': True, 'plan_sha256': plan_sha, 'outcomes_sha256': outcomes_sha,
                  'datasets': ['tenv'], 'tenv_test_memberships_previously_scored_in_step03': True, 'files': {}})
            destination = root / 'output/tenv'
            result = subprocess.run(
                [sys.executable, str(COLLECT), str(plan), plan_sha, str(outcomes),
                 '--publish-checkpoints', str(evaluation), str(destination)],
                capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            index = json.loads((destination / 'checkpoints.json').read_text())
            self.assertEqual([(record['identity']['fold'], record['identity']['seed']) for record in index['checkpoints']],
                             [(fold, seed) for fold in range(3) for seed in (17, 42, 2026)])
            for row in index['checkpoints']:
                original = root / row['directory']
                published = destination / row['directory']
                self.assertTrue((published / 'best.pt').samefile(original / 'best.pt'))
                self.assertEqual((published / 'replay.json').read_bytes(), (original / 'replay.json').read_bytes())

    def test_refuses_release_when_evaluation_is_incomplete(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = root / 'plan.json'
            plan_sha = write(plan, {'schema': 'step04-evaluation-plan-v1', 'stage': 'formal',
                                    'dataset': 'tenv', 'folds': [{}, {}, {}]})
            outcomes = root / 'outcomes.json'
            outcomes_sha = write(outcomes, {'schema': 'step04-evaluation-outcomes-v1',
                                             'plan_sha256': plan_sha, 'outcomes': []})
            evaluation = root / 'evaluation'
            write(evaluation / 'evaluation.json', {'schema': 'step04-evaluation-v1', 'stage': 'formal',
                  'complete': False, 'plan_sha256': plan_sha, 'outcomes_sha256': outcomes_sha,
                  'datasets': ['tenv'], 'tenv_test_memberships_previously_scored_in_step03': True})
            destination = root / 'output' / 'tenv'
            result = subprocess.run(
                [sys.executable, str(COLLECT), str(plan), plan_sha, str(outcomes),
                 '--publish-checkpoints', str(evaluation), str(destination)],
                capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('complete evaluated tenv matrix', result.stderr)
            self.assertFalse(destination.exists())


if __name__ == '__main__':
    unittest.main()
