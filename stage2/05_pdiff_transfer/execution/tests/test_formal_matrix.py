"""Fixture coverage for the fixed Step 05 plan and phase controls."""
from collections import Counter
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

import formal_matrix


class FormalMatrixTests(unittest.TestCase):
    def test_fixed_slots_and_worker_dependencies(self):
        slots, workers = formal_matrix.fixed_matrix()
        self.assertEqual(len({(row['regime'], row['fold'], row['arm'], row['seed']) for row in slots}), 36)
        self.assertEqual(Counter(row['operation'] for row in workers.values()),
                         {'fit': 21, 'source-fit': 3, 'predict': 21, 'tahoe-predict': 9})
        self.assertEqual(sum('reuse' in row for row in slots), 6)
        self.assertEqual(workers['cold-drug-tahoe-42-predict']['requires'], ['cold-drug-tahoe-42-fit'])
        self.assertEqual(workers['cold-drug-transfer-42-fit']['requires'], ['cold-drug-tahoe-42-predict'])
        self.assertEqual(workers['cold-drug-transfer-42-predict']['requires'], ['cold-drug-transfer-42-fit'])
        self.assertEqual(workers['warm-start-tahoe-17-predict']['requires'], [])

    def test_live_pause_then_explicit_resume_preserves_independent_seed_and_prediction_retry(self):
        from contextlib import nullcontext
        import subprocess
        import tempfile
        import threading
        import time
        from unittest.mock import patch

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            backend = formal_matrix.CappedScheduler(formal_matrix.transfer.REPO)
            scheduler = backend.scheduler
            slots, workers = formal_matrix.fixed_matrix()
            selected = {'cold-drug-tahoe-17-fit', 'cold-drug-tahoe-17-predict',
                        'cold-drug-tahoe-42-fit', 'cold-drug-tahoe-42-predict'}
            slots = [row for row in slots if row.get('prediction') in selected]
            workers = {key: row for key, row in workers.items() if key in selected}
            start_marker = root / 'started'
            retry_marker = root / 'retried'
            jobs = {
                'cold-drug-tahoe-17-fit': 'raise SystemExit(1)',
                'cold-drug-tahoe-42-fit': f'from pathlib import Path; import time; '
                                         f'Path({str(start_marker)!r}).touch(); time.sleep(.25)',
                'cold-drug-tahoe-42-predict': f'from pathlib import Path; p = Path({str(retry_marker)!r}); '
                                             'first = not p.exists(); p.touch(); raise SystemExit(int(first))',
            }
            for key, worker in workers.items():
                worker['queue_job'] = {'id': key, 'kind': 'cpu', 'cwd': temp, 'timeout_seconds': 2,
                                       'profile': {'host_peak_gib': .1, 'gpu_peak_gib': 0},
                                       'command': [sys.executable, '-c', jobs.get(key, 'raise SystemExit(99)')]}
            plan = {'slots': slots, 'workers': workers, 'attempt_limit': 2}
            observed, pause_requested = [], threading.Event()

            def pause_while_running():
                if not pause_requested.wait(timeout=4):
                    return
                observed.append(formal_matrix.control(root, 'pause', 'fixture relay'))

            def sample(record, pids=()):
                now = time.monotonic()
                return {'observable': True, 'available_gib': 40, 'gpu_free_gib': 15, 'gpu_used_gib': 1,
                        'gpu_utilization_percent': 0, 'disk_free_gib': 100, 'external_heavy': 0,
                        'external_peak_gib': 0, 'timestamp': time.time(), 'sample_started_monotonic': now,
                        'observed_monotonic': now}

            def launch(job, attempt):
                with Path(attempt['log']).open('w') as log:
                    process = subprocess.Popen(job['command'], cwd=temp, stdout=log, stderr=subprocess.STDOUT,
                                               start_new_session=True)
                if job['id'] == 'cold-drug-tahoe-42-fit':
                    pause_requested.set()
                return {'process': process, 'started': time.monotonic()}

            def verify(_plan, key, prior):
                if key == 'cold-drug-tahoe-42-predict':
                    self.assertEqual(prior['cold-drug-tahoe-42-fit']['status'], 'completed')
                return {'fixture': key}

            thread = threading.Thread(target=pause_while_running)
            thread.start()
            with (patch.object(scheduler, 'workstation_lock', nullcontext),
                  patch.object(scheduler, 'snapshot', sample),
                  patch.object(scheduler, 'launch', launch),
                  patch.object(scheduler, 'SAMPLE_SECONDS', .02)):
                self.assertFalse(formal_matrix.run_matrix(root, 'fixture', plan, backend, verify))
                thread.join(timeout=4)
                self.assertEqual(observed[0]['active'], ['cold-drug-tahoe-42-fit'])
                self.assertTrue(start_marker.exists())
                self.assertFalse((root / 'phases/tahoe-predictions/queue.json').exists())
                self.assertFalse(formal_matrix.run_matrix(root, 'fixture', plan, backend, verify))
                formal_matrix.control(root, 'resume', 'fixture explicit resume')
                self.assertTrue(formal_matrix.run_matrix(root, 'fixture', plan, backend, verify))

            inventory = formal_matrix.verify_inventory(root, 'fixture', plan)
            outcomes = {slot['seed']: slot for slot in inventory['slots']}
            self.assertEqual(outcomes[17]['status'], 'blocked')
            self.assertEqual(len(outcomes[17]['fit']['attempts']), 2)
            self.assertEqual(outcomes[17]['prediction']['failed_prerequisites'], ['cold-drug-tahoe-17-fit'])
            self.assertEqual(outcomes[42]['status'], 'completed')
            self.assertEqual(len(outcomes[42]['fit']['attempts']), 1)
            self.assertEqual([row['outcome'] for row in outcomes[42]['prediction']['attempts']],
                             ['failed', 'completed'])
            queue = formal_matrix.read_json(root / 'phases/tahoe-predictions/queue.json')
            self.assertEqual([row['id'] for row in queue['jobs']], ['cold-drug-tahoe-42-predict'])
            self.assertIn('fixture explicit resume', (root / 'control.jsonl').read_text())
            self.assertTrue((root / 'phases/tahoe-predictions/scheduler/admission.jsonl').is_file())
            path = root / 'prediction_inventory.json'
            inventory['slots'].pop()
            formal_matrix.write_json(path, inventory)
            with self.assertRaisesRegex(ValueError, 'inventory slots differ'):
                formal_matrix.verify_inventory(root, 'fixture', plan)
