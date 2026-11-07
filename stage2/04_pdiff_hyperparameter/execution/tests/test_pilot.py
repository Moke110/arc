"""Step 04 validation-only pilot artifact CLI."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

CLI = Path(__file__).resolve().parents[1] / 'script/pilot.py'


class PilotArtifactsTest(unittest.TestCase):
    def test_queue_pairs_cadences_for_every_tenv_fold(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            release = root / 'frozen-release'
            output = root / 'pilot'
            subprocess.run([sys.executable, str(CLI), 'profile', str(output), '--tenv-release', str(release)],
                           check=True, capture_output=True, text=True)
            queue = json.loads((output / 'profile' / 'queue.json').read_text())
            fit_jobs = [job for job in queue['jobs'] if job['kind'] == 'training']
            self.assertEqual(len(fit_jobs), 6)
            self.assertEqual(len([job for job in queue['jobs'] if job['kind'] == 'cpu']), 3)
            for fold in range(3):
                jobs = [job for job in fit_jobs if job['id'].startswith(f'fold-{fold}-')]
                recipes = [json.loads(Path(job['command'][-2]).read_text()) for job in jobs]
                self.assertEqual([recipe.pop('validation_interval') for recipe in recipes], [250, 500])
                self.assertEqual(recipes[0], recipes[1])
                self.assertEqual(recipes[0]['stage'], 'pilot')
                self.assertEqual(recipes[0]['fold'], fold)
                self.assertGreaterEqual(recipes[0]['max_steps'], 1000)
                self.assertGreater(recipes[0]['minimum_steps'], recipes[0]['max_steps'])

    def test_cadence_requires_matched_traces_and_uses_measured_overhead(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'pilot'
            self.make_profile(output)
            report = subprocess.run([sys.executable, str(CLI), 'calibrate', str(output)],
                                    capture_output=True, text=True)
            self.assertEqual(report.returncode, 0, report.stderr)
            selection = json.loads((output / 'cadence.json').read_text())
            self.assertEqual([fold['validation_interval'] for fold in selection['folds']], [500, 250, 250])
            self.assertEqual(selection['folds'][0]['profiles'][0]['overhead_fraction'], 0.32)

    def test_divergent_matched_step_rejects_calibration(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'pilot'
            self.make_profile(output)
            trace = output / 'profile/runs/fold-0-v250/run.json'
            run = json.loads(trace.read_text())
            event = run['history'][1]
            event['validation_loss'] -= 0.01
            event['batch_losses'] = [event['validation_loss']]
            trace.write_text(json.dumps(run))
            report = subprocess.run([sys.executable, str(CLI), 'calibrate', str(output)],
                                    capture_output=True, text=True)
            self.assertNotEqual(report.returncode, 0)
            self.assertIn('pilot validation trajectories diverge at 500', report.stderr)
            self.assertFalse((output / 'calibration').exists())

    def test_freezes_recipes_from_plateau_and_reports_censored_trajectory(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'pilot'
            self.make_profile(output)
            subprocess.run([sys.executable, str(CLI), 'calibrate', str(output)],
                           check=True, capture_output=True, text=True)
            self.make_phase(output / 'calibration', {
                'fold-0-v500': lambda step: 10 - step / 2000,
                'fold-1-v250': lambda step: 10 - min(step, 4000) / 2000,
                'fold-2-v250': lambda step: 10 - min(step, 6000) / 2000,
            })
            report = subprocess.run([sys.executable, str(CLI), 'select', str(output)],
                                    capture_output=True, text=True)
            self.assertEqual(report.returncode, 0, report.stderr)
            frozen = json.loads((output / 'selection' / 'selection.json').read_text())
            fold0, fold1, fold2 = frozen['folds']
            self.assertEqual(fold0['calibration_status'], 'right_censored')
            self.assertEqual(fold1['selected_candidate']['stop_step'], 5000)
            self.assertEqual(fold2['selected_candidate']['stop_step'], 7000)
            for fold in (fold0, fold1, fold2):
                self.assertEqual(set(fold['formal_recipes']), {'17', '42', '2026'})
                for seed, entry in fold['formal_recipes'].items():
                    recipe = json.loads((output / 'selection' / entry['path']).read_text())
                    self.assertEqual(recipe['stage'], 'formal')
                    self.assertEqual(recipe['seed'], int(seed))
                    self.assertEqual(recipe['min_delta'], fold['selected_candidate']['min_delta'])
            compact = Path(temporary) / 'compact'
            (compact / 'calibration').mkdir(parents=True)
            shutil.copy2(output / 'cadence.json', compact)
            shutil.copy2(output / 'calibration/plan.json', compact / 'calibration')
            shutil.copytree(output / 'selection', compact / 'selection')
            verify = [sys.executable, str(CLI), 'verify', str(compact)]
            subprocess.run(verify, check=True, capture_output=True, text=True)
            recipe = compact / 'selection/fold-0/seed-17.json'
            edited = json.loads(recipe.read_text())
            edited['minimum_steps'] -= 1
            recipe.write_text(json.dumps(edited))
            rejected = subprocess.run(verify, capture_output=True, text=True)
            self.assertIn('formal recipes do not match', rejected.stderr)

    def make_profile(self, output):
        requests = []
        for fold in range(3):
            request = output.parent / f'fold-{fold}.json'
            cases = request.with_suffix('.csv')
            with cases.open('w') as stream:
                writer = csv.writer(stream)
                writer.writerows([['group_id'], ['one']])
            request.write_text(json.dumps({'identity': {'dataset': 'tenv', 'fold': fold},
                                           'predictor_files': {'validation_cases': str(cases)}}))
            requests.append(str(request))
        subprocess.run([sys.executable, str(CLI), 'profile', str(output), '--requests', *requests],
                       check=True, capture_output=True, text=True)
        self.make_phase(output / 'profile', {
            f'fold-{fold}-v{interval}': (lambda step: 10 - step / 2000)
            for fold in range(3) for interval in (250, 500)
        })

    def make_phase(self, phase, monitors):
        queue = json.loads((phase / 'queue.json').read_text())
        status = {}
        for job in queue['jobs']:
            job_id = job['id']
            command = job['command']
            recipe = json.loads(Path(command[-2]).read_text())
            path = Path(command[-1]) / 'run.json'
            path.parent.mkdir(parents=True)
            request = json.loads(Path(command[-3]).read_text())
            interval = recipe['validation_interval']
            overhead = 4 if recipe['fold'] == 0 else 1
            history = []
            best_step, best_loss = None, float('inf')
            for step in range(interval, recipe['max_steps'] + 1, interval):
                loss = monitors[job_id](step)
                refreshed = loss < best_loss
                if refreshed:
                    best_step, best_loss = step, loss
                history.append({'step': step, 'validation_loss': loss, 'batch_losses': [loss],
                                'batch_sizes': [1], 'validation_seconds': overhead,
                                'checkpoint_seconds': 0, 'checkpoint_refreshed': refreshed,
                                'training_seconds': 100 * step / recipe['max_steps'], 'training_loss': 5})
            digest = lambda data: hashlib.sha256(
                json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            run = {'identity': {**request['identity'], 'stage': 'pilot', 'seed': 17},
                   'status': 'completed', 'stop_reason': 'max_steps', 'optimizer_steps': recipe['max_steps'],
                   'recipe': {**recipe, 'batch_groups': 8, 'population_size': 3},
                   'monitor': 'validation_validation_loss_epoch', 'history': history,
                   'best_step': best_step, 'best_validation_loss': best_loss, 'training_seconds': 100,
                   'wall_seconds': 160, 'reload_seconds': 2, 'upstream_commit': 'fixture',
                   'reload': {'validation_loss': best_loss, 'prediction_max_abs_error': 0},
                   'resources': {'peak_rss_kib': 1, 'peak_gpu_reserved_bytes': 1},
                   'checkpoint_sha256': 'a' * 64, 'model_state_sha256': 'b' * 64,
                   'environment': {}, 'validation_batching': {'group_ids': ['one'],
                   'group_ids_sha256': digest(['one']), 'batch_groups': 8, 'population_size': 3,
                   'random_seed': 100017, 'remainder': 'keep and weight as one batch'}}
            run['recipe_sha256'] = digest(run['recipe'])
            path.write_text(json.dumps(run))
            status[job_id] = {'status': 'completed', 'identity': job, 'attempts': [{
                'outcome': 'completed', 'returncode': 0, 'peak_host_rss_gib': 1, 'peak_workstation_gpu_used_gib': 1,
                'started_at': 0, 'finished_at': 160, 'memory_max_gib': 10, 'memory_swap_max_gib': 8
            }]}
        path = phase / 'scheduler' / 'status.json'
        path.parent.mkdir()
        path.write_text(json.dumps({'jobs': status}))


if __name__ == '__main__':
    unittest.main()
