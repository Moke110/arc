"""Historical metadata remains in Git; frozen predictor proof stays authenticated."""
import copy
import hashlib
import json
from pathlib import Path
import unittest
from unittest.mock import patch
from stage2.evaluation import compatibility
from stage2.evaluation.history import REPO, original_json, original_bytes


class HistoryTest(unittest.TestCase):
    def test_original_publication_proves_existing_checkpoint_index_without_restoring_scores(self):
        path = REPO / 'stage2/05_pdiff_transfer/execution/formal_record/publication.json'
        self.assertFalse(path.exists())
        previous = original_json(path)
        checkpoints = path.with_name('checkpoints.json')
        self.assertEqual(previous['files']['checkpoints.json'], hashlib.sha256(checkpoints.read_bytes()).hexdigest())
        with self.assertRaisesRegex(ValueError, 'digest differs'):
            original_bytes(path, '0' * 64)
        self.assertFalse(path.exists())

    def test_only_reviewed_runtime_and_original_frozen_run_can_bridge_source_identity(self):
        inventory = json.loads((REPO / 'stage2/evaluation/inventory.json').read_text())
        item = next(item for item in inventory['runs'] if item['step'] == 3)
        path = REPO / item['run']['path']
        previous = json.loads(path.read_text())['identity']
        current = {**previous, 'source_sha256': 'new-reporting-source'}
        self.assertTrue(compatibility.accepts_frozen_run(path, previous, current))
        altered = {**current, 'input_sha256': 'changed-input'}
        self.assertFalse(compatibility.accepts_frozen_run(path, previous, altered))
        self.assertFalse(compatibility.accepts_frozen_run(path.with_name('unknown-run.json'), previous, current))
        with patch.object(compatibility.hashlib, 'sha256') as digest:
            digest.return_value.hexdigest.return_value = 'changed-runtime'
            self.assertFalse(compatibility.accepts_frozen_run(path, previous, current))

    def test_original_completion_barrier_uses_authenticated_identity_bridge(self):
        import importlib.util
        import sys
        script = REPO / 'stage2/03_modes_tenv/execution/script/benchmark.py'
        spec = importlib.util.spec_from_file_location('frozen_benchmark', script)
        module = importlib.util.module_from_spec(spec)
        with patch.object(sys, 'path', [str(script.parent), *sys.path]):
            spec.loader.exec_module(module)
        output = REPO / 'stage2/03_modes_tenv/output/formal'
        def identity(inputs, model, fold, seed, command):
            path = module.run_directory(output, model, fold, seed) / 'run.json'
            previous = json.loads(path.read_text())['identity']
            return {**previous, 'source_sha256': 'new-reporting-source'}
        with patch.object(module, 'run_identity', identity), patch.object(module, 'validate_run') as validate, \
             patch.object(module, 'write_json'):
            result = module.freeze_predictions(REPO / 'stage2/03_modes_tenv/output/input', output, ['no_change'])
        self.assertTrue(result['sealed'])
        self.assertEqual(validate.call_count, 9)


if __name__ == '__main__':
    unittest.main()
