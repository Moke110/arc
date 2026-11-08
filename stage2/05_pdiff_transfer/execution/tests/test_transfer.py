"""Checks at the Transfer artifact command boundary."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


CLI = Path(__file__).resolve().parents[1] / 'script/transfer.py'


class TransferCLI(unittest.TestCase):
    def test_verification_rejects_an_unsealed_plan_before_accessing_inputs(self):
        result = subprocess.run(
            [sys.executable, str(CLI), 'verify', str(CLI.parent), '--plan-sha256', '0' * 64],
            text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn('Transfer verify rejected', result.stderr)


class SealedTransferFit(unittest.TestCase):
    def setUp(self):
        self.addCleanup(sys.path.remove, str(CLI.parent))
        sys.path.insert(0, str(CLI.parent))
        from transfer import verify_transfer_fit
        from artifacts import sha256, write_json
        from training_inputs import json_hash

        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.slot = Path(temporary.name)
        self.fit = self.slot / 'fit'
        self.fit.mkdir()
        self.verify_transfer_fit = verify_transfer_fit
        self.write_json, self.sha256 = write_json, sha256
        self.plan_sha = 'f' * 64
        self.identity = {'dataset': 'tenv', 'arm': 'Transfer', 'regime': 'cold-cm-drug',
                         'fold': 0, 'stage': 'formal', 'seed': 17}
        recipe = {'fold': 0, 'stage': 'formal', 'seed': 17}
        self.item = {'config': recipe}
        self.lineage = {'bridge_sha256': 'b' * 64}
        (self.fit / 'best.pt').write_bytes(b'selected EMA')
        self.run = {
            'status': 'completed', 'monitor': 'validation_validation_loss_epoch',
            'identity': self.identity, 'recipe': recipe, 'recipe_sha256': json_hash(recipe),
            'history': [{'step': 12, 'validation_loss': 1.0}], 'best_step': 12,
            'best_validation_loss': 1.0, 'reload': {'validation_loss': 1.0, 'prediction_max_abs_error': 0},
            'initialization': self.lineage, 'checkpoint_sha256': sha256(self.fit / 'best.pt'),
            'training_initial_state': {
                'optimizer_state_entries': 0, 'scheduler_start_step': 0, 'ema_start_step': 0,
                'mapped_dose_rows_trainable': True, 'mapped_dose_gradient_observed': True,
            },
        }
        self.seal_run()

    def seal_run(self):
        self.write_json(self.fit / 'run.json', self.run)
        self.write_json(self.fit / 'fit.json', {
            'plan_sha256': self.plan_sha, 'run_sha256': self.sha256(self.fit / 'run.json'),
        })

    def verify(self):
        return self.verify_transfer_fit(self.slot, self.plan_sha, self.identity, self.item, self.lineage)

    def test_resealed_fit_cannot_fill_another_seed_slot(self):
        self.assertEqual(self.verify()['identity']['seed'], 17)
        self.run['identity'] = {**self.identity, 'seed': 42}
        self.seal_run()
        with self.assertRaisesRegex(ValueError, 'formal identity or recipe'):
            self.verify()

    def test_substituted_checkpoint_cannot_fill_the_slot(self):
        self.assertEqual(self.verify()['identity']['seed'], 17)
        (self.fit / 'best.pt').write_bytes(b'other checkpoint')
        with self.assertRaisesRegex(ValueError, 'frozen artifact path/hash mismatch: best.pt'):
            self.verify()
