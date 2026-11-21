"""Behavior at the successor plan and selection seams."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / 'output/api'))
from tsg_round1.inputs import sha256
from tsg_round1.matrix import (
    check_registration, lock_selection, plan, read_json, register, run_path, select,
    validate_origins,
)
from tenv5k.attempt_identity import attempt_binding


class MatrixTests(unittest.TestCase):
    def test_formal_plan(self):
        origins = plan()
        self.assertEqual(len(origins), 84)
        self.assertEqual(len({tuple(origin.values()) for origin in origins}), 84)
        self.assertEqual({origin['seed'] for origin in origins}, {42})
        self.assertEqual(len(plan('cold-drug')), 42)
        self.assertEqual(len(plan('cold-drug-dose')), 42)
        with self.assertRaises(ValueError):
            plan('final')

    def test_cold_drug_selects_by_validation_then_size(self):
        origins = [origin for origin in plan('cold-drug')
                   if origin['source'] == 'tenv' and origin['fold'] == 0]
        scores = {origin['method']: (0.2, 100) for origin in origins}
        scores['film'] = (0.5, 800)
        scores['scaling'] = (0.5, 100)
        self.assertEqual(select(scores, 'cold-drug'), 'scaling')
        self.assertIsNone(select(scores, 'cold-drug-dose'))
        with self.assertRaises(ValueError):
            select({**scores, 'extra': (1., 1)}, 'cold-drug')

    def test_formal_acceptance_rejects_missing_duplicate_and_foreign_origins(self):
        expected = plan('cold-drug')
        validate_origins(expected, list(reversed(expected)))
        for found in (expected[:-1], expected[:-1] + expected[:1], expected[:-1] + plan('cold-drug-dose')[:1]):
            with self.assertRaisesRegex(ValueError, 'formal origins'):
                validate_origins(expected, found)

    def test_registration_binds_source_generation_and_origin(self):
        origin = plan()[0]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'input' / origin['source']
            source.mkdir(parents=True)
            manifest = source / 'manifest.json'
            manifest.write_text('accepted')
            registration_path = register(root / 'matrix', source.parent, origin, 'cuda')
            registration = read_json(registration_path)
            check_registration(registration, source, origin, 'cuda')
            with self.assertRaises(FileExistsError):
                register(root / 'matrix', source.parent, origin, 'cuda')
            with self.assertRaisesRegex(ValueError, 'registration identity'):
                check_registration(registration, source, dict(origin, seed=17), 'cuda')
            manifest.write_text('different generation')
            with self.assertRaisesRegex(ValueError, 'registration identity'):
                check_registration(registration, source, origin, 'cuda')

    def test_selection_locks_only_complete_validation_evidence(self):
        origins = [origin for origin in plan('cold-drug')
                   if origin['source'] == 'tenv' and origin['fold'] == 0]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = root / 'input' / 'tenv'
            inputs.mkdir(parents=True)
            (inputs / 'manifest.json').write_text('accepted')
            for origin in origins:
                run = run_path(root, origin)
                run.mkdir(parents=True)
                registered = read_json(register(root, inputs.parent, origin, 'cuda'))
                bound = attempt_binding(registered)
                scores = {'pcc': .75 if origin['method'] == 'film' else .5}
                terminal = dict(status='validation-complete', smoke=False, parameter_count=100,
                                training_process_token='training',
                                validation_replay=dict(process_token='validation',
                                                       attempt_nonce=registered['attempt_nonce']))
                for name, value in [('config.json', dict(origin, device='cuda', execution_binding=bound)),
                                    ('registration.json', registered),
                                    ('validation_metrics.json', scores),
                                    ('terminal.json', terminal)]:
                    (run / name).write_text(json.dumps(value))
                names = ['checkpoint.pt', 'environment.json', 'learning_curve.csv',
                         'validation_identity.csv', 'validation_per_sample.csv', 'validation_prediction.npy']
                for name in names:
                    (run / name).write_text('validation payload')
                evidence = dict(schema='tsg-dose-round1-axis-local-validation-v1',
                                input_manifest_sha256=sha256(inputs / 'manifest.json'),
                                parameter_count=100, outputs={path.name: sha256(path) for path in run.iterdir()
                                                              if path.name != 'terminal.json'})
                (run / 'validation_manifest.json').write_text(json.dumps(evidence))
            missing = run_path(root, origins[-1]) / 'validation_manifest.json'
            saved = missing.read_bytes()
            missing.unlink()
            with self.assertRaises(FileNotFoundError):
                lock_selection(root, inputs, origins[0])
            missing.write_bytes(saved)
            terminal_path = run_path(root, origins[-1]) / 'terminal.json'
            terminal_bytes = terminal_path.read_bytes()
            completed = dict(read_json(terminal_path), status='complete')
            terminal_path.write_text(json.dumps(completed))
            with self.assertRaisesRegex(ValueError, 'lifecycle mismatch'):
                lock_selection(root, inputs, origins[0])
            terminal_path.write_bytes(terminal_bytes)
            locked = lock_selection(root, inputs, origins[0])
            self.assertEqual(read_json(locked)['selected_method'], 'film')
            self.assertEqual(lock_selection(root, inputs, origins[0]), locked)
            (run_path(root, origins[-1]) / 'validation_metrics.json').write_text('{"pcc":0.9}')
            with self.assertRaisesRegex(ValueError, 'validation artifact changed'):
                lock_selection(root, inputs, origins[0])
