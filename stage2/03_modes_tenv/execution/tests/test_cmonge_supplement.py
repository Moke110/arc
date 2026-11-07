"""Compact supplemental-account and phase-watchdog tests without model or data."""
import json
import os
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest.mock import Mock, patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
import cmonge_supplement as runtime


class PhaseTests(unittest.TestCase):
    def test_valid_transition_extends_only_inference(self):
        clock = runtime.PhaseDeadline(50)
        self.assertIsNone(clock.observe(100, dict(phase='training', training_seconds=49)))
        self.assertIsNone(clock.observe(14449, dict(phase='inference', training_seconds=14300)))
        self.assertEqual(clock.deadline, 16249)
        self.assertIsNone(clock.observe(16000, dict(phase='inference', training_seconds=14300)))
        self.assertEqual(clock.observe(16249, None), 'inference_timeout')
        self.assertEqual(clock.timing(16249)['observed_inference_seconds'], 1800)

    def test_late_transition_cannot_authorize_extra_wall(self):
        clock = runtime.PhaseDeadline(50)
        self.assertEqual(clock.observe(14450, dict(phase='inference', training_seconds=14300)),
                         'training_timeout')
        self.assertIsNone(clock.transition)

    def test_invalid_progress_refuses_credit(self):
        for value in (True, float('inf'), float('nan'), -1, 100):
            with self.subTest(value=value), self.assertRaises(ValueError):
                runtime.PhaseDeadline(0).observe(10, dict(phase='inference', training_seconds=value))
        clock = runtime.PhaseDeadline(0)
        clock.observe(6, dict(phase='training', training_seconds=5))
        for progress in (None, dict(phase='inference', training_seconds=4)):
            with self.assertRaises(ValueError):
                clock.observe(7, progress)
        clock.observe(8, dict(phase='inference', training_seconds=5))
        with self.assertRaises(ValueError):
            clock.observe(9, dict(phase='training', training_seconds=6))

    def test_zero_exit_requires_valid_transition_and_deadline(self):
        with tempfile.TemporaryDirectory() as temp:
            progress = Path(temp) / 'progress.json'
            transition = Path(temp) / 'transition.json'
            worker = Mock()
            worker.poll.return_value = 0
            stop = threading.Event()
            with patch.object(runtime.time, 'monotonic', return_value=10):
                self.assertEqual(runtime.watch(worker, runtime.PhaseDeadline(0), progress, stop, transition)[:2],
                                 ('failed', 'worker_exit_without_inference'))
                worker.poll.return_value = 124
                self.assertEqual(runtime.watch(worker, runtime.PhaseDeadline(0), progress, stop, transition)[:2],
                                 ('timeout', 'native_training_timeout'))
                worker.poll.return_value = 0
                runtime.write(progress, dict(phase='inference', training_seconds=3))
                self.assertEqual(runtime.watch(worker, runtime.PhaseDeadline(0), progress, stop, transition)[:2],
                                 ('completed', 'worker_zero_exit'))
            with patch.object(runtime.time, 'monotonic', return_value=runtime.TRAIN):
                self.assertEqual(runtime.watch(worker, runtime.PhaseDeadline(0), progress, stop, transition)[:2],
                                 ('timeout', 'training_timeout'))


class AccountTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.output = self.root / 'output'
        self.record = self.root / 'record'
        self.lock = self.root / '.scheduler.lock'
        self.lock.touch()
        patches = patch.multiple(runtime, OUTPUT=self.output, RECORD=self.record, LOCK=self.lock)
        patches.start()
        self.addCleanup(patches.stop)
        self.directory = runtime.run_directory(self.output, *runtime.SLOT)
        self.account = self.output / 'account'
        self.policy = dict(identity={'model': 'cmonge', 'fold': 0, 'seed': 17},
                           caps=dict(training_seconds=14400, inference_seconds=1800),
                           job={'command': ['/env/python', '/frozen/cmonge.py']},
                           original={'status': 'timeout'}, account=str(self.account),
                           terminal=str(self.record / 'terminal.json'))
        self.cases = self.root / 'cases.csv'
        self.cases.write_text('case_id\nheld-out-case\n')
        self.request = dict(identity=self.policy['identity'], output_directory=str(self.directory),
                            predictor_files={'test_cases': str(self.cases)})
        self.digest = 'a' * 64

    def reserve(self):
        fd = runtime.lock_fd()
        self.addCleanup(os.close, fd)
        return runtime.reserve(self.policy, self.digest, self.request, fd, {}, {})

    def test_canonical_lock_and_permanent_one_shot(self):
        state = self.reserve()
        self.assertEqual(state['reservation_count'], 1)
        self.assertEqual(state['credit_available_seconds'], 0)
        self.assertIn('MemoryMax=10G', state['scope_command'])
        self.assertIn('MemorySwapMax=8G', state['scope_command'])
        self.assertIn('RuntimeMaxSec=16200', state['scope_command'])
        self.assertEqual(state['executed_argv'][-2:], ['--request', str(self.account / 'request-0.json')])
        self.assertEqual(runtime.read_json(self.account / 'request-0.json'),
                         dict(self.request, prior_training_seconds=0))
        with self.assertRaises(BlockingIOError):
            runtime.lock_fd()
        with self.assertRaisesRegex(ValueError, 'already claimed'):
            runtime.untouched()

    def test_symlink_and_hardlink_refused(self):
        self.output.symlink_to(self.root / 'absent')
        with self.assertRaisesRegex(ValueError, 'already claimed'):
            runtime.untouched()
        extra = self.root / 'hardlink'
        os.link(self.lock, extra)
        with self.assertRaisesRegex(ValueError, 'single-link'):
            runtime.lock_fd()

    def test_scope_cleanup_requires_inactive_proof(self):
        worker = Mock()
        worker.poll.return_value = 0
        with patch.object(runtime, 'systemctl') as systemctl, patch.object(runtime, 'unit_properties') as properties:
            systemctl.return_value.returncode = 0
            properties.return_value = dict(ActiveState='active', LoadState='loaded')
            with self.assertRaisesRegex(ValueError, 'scope exit unverified'):
                runtime.cleanup(worker, 'test.scope', None)
            properties.return_value = dict(ActiveState='inactive', LoadState='not-found')
            self.assertEqual(runtime.cleanup(worker, 'test.scope', None)['active_state'], 'inactive')

    def complete_checkpoint(self):
        events = [dict(iteration=i, validation_loss=1., validation_conditions_evaluated=138)
                  for i in range(50, 1001, 50)]
        history = dict(identity=self.policy['identity'], protocol=runtime.PROTOCOL,
                       snapshot_semantics='prediction-only-not-resumable', history=events)
        runtime.write(self.directory / 'history.json', history)
        for name in ('best.msgpack', 'latest.msgpack'):
            (self.directory / name).write_bytes(b'fixture')
        runtime.write(self.directory / 'progress.json', dict(phase='inference', training_seconds=15))
        runtime.write(self.directory / 'case_ids.json', ['held-out-case'])
        np.save(self.directory / 'prediction.npy', np.zeros((1, 18164), dtype=np.float32))
        return dict(history, prior_training_seconds=0, gene_count=18164, backend='gpu',
                    validation_groups=138, batch_size=512, maximum_iterations=10000,
                    qualification_only=False, training_seconds=15, inference_seconds=1,
                    blobs={name: runtime.sha256(self.directory / name) for name in
                           ('best.msgpack', 'latest.msgpack', 'history.json')},
                    early_stopping=dict(minimum_iterations=500, patience=10, relative_min_delta=.001,
                                        earliest_stop_iteration=1000),
                    iterations=1000, best_iteration=50, best_validation_loss=1.,
                    stop_reason='validation_patience')

    def test_completed_proof_and_tampering(self):
        state = self.reserve()
        checkpoint = self.complete_checkpoint()
        runtime.write(self.directory / 'checkpoint.json', checkpoint)
        state.update(status='completed', stop_reason='worker_zero_exit', returncode=0,
                     timing=dict(adapter_training_seconds=15, worker_wall_seconds=20,
                                 pre_inference_watchdog_seconds=18, observed_inference_seconds=2),
                     cleanup=dict(scope=state['scope'], active_state='inactive'))
        terminal = runtime.publish_terminal(self.policy, self.digest, state, self.request)
        self.assertEqual(terminal['status'], 'completed')
        with patch.object(runtime, 'checked_policy', return_value=self.policy), \
             patch.object(runtime, 'request_for', return_value=self.request):
            terminal_sha = runtime.sha256(self.record / 'terminal.json')
            args = (None, None, None, self.digest, terminal_sha)
            self.assertEqual(runtime.verify_terminal(*args)[1], terminal)
            (self.directory / 'best.msgpack').write_bytes(b'different')
            with self.assertRaises(ValueError):
                runtime.verify_terminal(*args)

    def test_timeout_keeps_partial_checkpoint_out_of_scoring(self):
        state = self.reserve()
        state.update(status='timeout', stop_reason='training_timeout', returncode=-9,
                     timing=dict(adapter_training_seconds=14390, worker_wall_seconds=14400,
                                 pre_inference_watchdog_seconds=14400, observed_inference_seconds=None),
                     cleanup=dict(scope=state['scope'], active_state='inactive'))
        (self.directory / 'best.msgpack').write_bytes(b'partial')
        terminal = runtime.publish_terminal(self.policy, self.digest, state, self.request)
        self.assertEqual(terminal['status'], 'timeout')
        self.assertIn('best.msgpack', terminal['files'])
        self.assertNotIn('artifacts', runtime.read_json(self.directory / 'run.json'))

    def test_bad_validation_count_rejected(self):
        self.reserve()
        checkpoint = self.complete_checkpoint()
        checkpoint['history'][0]['validation_conditions_evaluated'] = 137
        runtime.write(self.directory / 'history.json',
                      {k: checkpoint[k] for k in ('protocol', 'identity', 'history', 'snapshot_semantics')})
        checkpoint['blobs']['history.json'] = runtime.sha256(self.directory / 'history.json')
        runtime.write(self.directory / 'checkpoint.json', checkpoint)
        with self.assertRaisesRegex(ValueError, 'excluded conditions'):
            runtime.validate_completion(self.directory, self.request,
                dict(timing=dict(adapter_training_seconds=15, worker_wall_seconds=20)))


if __name__ == '__main__':
    unittest.main()
