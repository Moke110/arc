"""Formal launcher boundaries stop on failed origins and never retrain completed phases."""
import importlib.util
from pathlib import Path
import unittest

path = Path(__file__).resolve().parents[1] / 'axis_local/formal.py'
spec = importlib.util.spec_from_file_location('axis_local_formal', path)
formal = importlib.util.module_from_spec(spec)
spec.loader.exec_module(formal)


class PhaseTests(unittest.TestCase):
    def test_resume_runs_only_missing_phase(self):
        self.assertEqual(formal.next_action(None, 'validation'), 'train')
        self.assertEqual(formal.next_action('validation-trained-awaiting-replay', 'validation'), 'replay')
        self.assertIsNone(formal.next_action('validation-complete', 'validation'))
        self.assertEqual(formal.next_action('validation-complete', 'test'), 'evaluate')
        self.assertEqual(formal.next_action('test-awaiting-replay', 'test'), 'replay')
        self.assertIsNone(formal.next_action('complete', 'test'))
        for status in ('failed', 'interrupted'):
            with self.assertRaisesRegex(ValueError, status):
                formal.next_action(status, 'validation')
