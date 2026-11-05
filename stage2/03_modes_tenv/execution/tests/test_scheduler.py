"""Resource policy is tested with telemetry, never by exhausting the host."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'script' / 'schedule.py'


class SchedulerContract(unittest.TestCase):
    def policy(self, snapshot, profile=None, active=None, max_active_jobs=2):
        with tempfile.TemporaryDirectory() as td:
            source = Path(td) / 'request.json'
            source.write_text(json.dumps(dict(snapshot=snapshot, profile=profile,
                                             active=active or [], max_active_jobs=max_active_jobs)))
            p = subprocess.run([sys.executable, str(SCRIPT), 'admit', str(source)],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            return json.loads(p.stdout)

    def test_unknown_profile_starts_alone_and_respects_external_work(self):
        snapshot = dict(available_gib=30, gpu_free_gib=12, external_heavy=0,
                        external_peak_gib=0, observable=True)
        self.assertTrue(self.policy(snapshot)['admit'])
        snapshot['external_heavy'] = 1
        snapshot['external_peak_gib'] = 3
        self.assertFalse(self.policy(snapshot)['admit'])

    def test_known_profiles_need_headroom_and_host_visibility(self):
        s = dict(available_gib=24, gpu_free_gib=10, external_heavy=1,
                 external_peak_gib=8, observable=True)
        self.assertTrue(self.policy(s, dict(host_peak_gib=8, gpu_peak_gib=4))['admit'])
        self.assertFalse(self.policy(s, dict(host_peak_gib=2),
                         active=[dict(profile=dict(host_peak_gib=2))], max_active_jobs=1)['admit'])
        self.assertFalse(self.policy(s, dict(host_peak_gib=12, gpu_peak_gib=4))['admit'])
        s['observable'] = False
        self.assertFalse(self.policy(s, dict(host_peak_gib=1))['admit'])

    def test_cumulative_budget_and_independent_retry_limits(self):
        with tempfile.TemporaryDirectory() as td:
            state = Path(td) / 'state.json'
            def transition(reason, seconds):
                p = subprocess.run([sys.executable, str(SCRIPT), 'transition',
                                    str(state), '--reason', reason,
                                    '--seconds', str(seconds)], capture_output=True, text=True)
                self.assertEqual(p.returncode, 0, p.stderr)
                return json.loads(p.stdout)
            first = transition('environment', 20)
            self.assertEqual(first['status'], 'pending')
            self.assertEqual(first['environment_retries'], 1)
            second = transition('oom', 30)
            self.assertEqual(second['batch_scale'], 0.5)
            self.assertEqual(second['training_seconds'], 50)
            self.assertEqual(transition('oom', 10)['status'], 'failed')
            # A terminal state cannot be reopened by another launch.
            self.assertEqual(transition('environment', 0)['status'], 'failed')
            state.unlink()
            self.assertEqual(transition('oom', 14400)['status'], 'timeout')
            state.unlink()
            inference = transition('inference_timeout', 0)
            self.assertEqual(inference['status'], 'timeout')
            self.assertEqual(inference['training_seconds'], 0)

    def test_recovered_terminal_states_satisfy_prediction_barrier(self):
        # Durable scheduler state is the recovery input after a dead worker or
        # scheduler shutdown. Both must become terminal artifact outcomes.
        from test_pipeline import fixture
        from artifacts import run_identity, write_json, run_directory
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); fixture(root/'input')
            jobs = []
            for fold in range(3):
                for seed in (17,42,2026):
                    identity = run_identity(root/'input','no_change',fold,seed)
                    key = f'no_change_fold{fold}_seed{seed}'
                    dead = fold == 0 and seed == 17
                    state = dict(identity=identity,status='running' if dead else 'failed',
                                 pid=2147483647,training_seconds=12.,stop_reason='scheduler_interrupted')
                    write_json(root/'output/scheduling'/f'{key}.json',state)
                    write_json(run_directory(root/'output','no_change',fold,seed)/'run.json',
                               dict(identity=identity,status='running'))
                    jobs.append(dict(model='no_change',fold=fold,seed=seed))
            write_json(root/'queue.json',dict(input_root=str(root/'input'),output_root=str(root/'output'),jobs=jobs))
            # Isolate the installed CLI's workstation lock from real experiments.
            import shutil
            isolated = root/'isolated/script'; isolated.mkdir(parents=True)
            for module in SCRIPT.parent.glob('*.py'):
                shutil.copyfile(module,isolated/module.name)
            result = subprocess.run([sys.executable,str(isolated/'schedule.py'),'run',str(root/'queue.json')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            benchmark = SCRIPT.with_name('benchmark.py')
            result = subprocess.run([sys.executable,str(benchmark),'seal','--input',str(root/'input'),
                                     '--output',str(root/'output'),'--models','no_change'],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            manifest = json.loads((root/'output/predictions.json').read_text())
            self.assertTrue(all(run['status']=='failed' for run in manifest['runs']))

    def test_pressure_stops_only_newest_owned_task(self):
        with tempfile.TemporaryDirectory() as td:
            source = Path(td) / 'pressure.json'
            source.write_text(json.dumps(dict(
                history=[dict(available_gib=11, gpu_free_gib=2),
                         dict(available_gib=10, gpu_free_gib=1)],
                tasks=[dict(id='older', started=1, owned=True),
                       dict(id='unrelated', started=3, owned=False),
                       dict(id='newer', started=2, owned=True)])))
            p = subprocess.run([sys.executable, str(SCRIPT), 'pressure', str(source)],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(json.loads(p.stdout)['stop'], 'newer')


if __name__ == '__main__':
    unittest.main()
