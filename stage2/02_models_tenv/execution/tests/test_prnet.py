"""PRnet at the agreed isolated fit/save/reload/predict artifact seam (#163).

The test drives the real command in the pinned PRnet environment on a small
synthetic release of the documented input contract.  It asserts the artifact
boundary rather than internals: a completed run persists a checkpoint, a
reloaded artifact reproduces predictions exactly, every case identity keeps its
order and shape, the early-stopping record names the ADR 0063 rules, and an
over-budget run records an honest terminal status instead of a partial
prediction.
"""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

STEP = Path(__file__).resolve().parents[2]
REPO = STEP.parents[1]
ADAPTER = STEP / 'execution/script/prnet_fit_predict.py'
BASE_CONFIG = STEP / 'execution/models/prnet/config.json'
PYTHON = '/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python'
N_GENES = 18164
MORGAN_BITS = 2048
CASE_FIELDS = ['case_id', 'drug_key', 'dose_um', 'group_id', 'control_group',
               'source_row', 'morgan_index']


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path, records):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def write_release(root, rng, *, n_train=6, n_validation=2, n_test=2,
                  delta_value=None):
    release = root / 'release'
    release.mkdir()
    drugs = [{'drug_key': f'drug{i}', 'morgan_index': i} for i in range(3)]
    write_csv(release / 'drugs.csv', drugs)
    np.save(release / 'morgan.npy',
            rng.integers(0, 2, (len(drugs), MORGAN_BITS), dtype='uint8'))
    files = {'drugs': 'drugs.csv', 'morgan': 'morgan.npy'}
    roles = [('train', 0, n_train), ('validation', 1, n_validation), ('test', 2, n_test)]
    for role, drug, count in roles:
        cases = [{'case_id': f'{role}{index}', 'drug_key': f'drug{drug}',
                  'dose_um': float(1 + index), 'group_id': f'group-{role}-{index}',
                  'control_group': 'synthetic', 'source_row': index,
                  'morgan_index': drug} for index in range(count)]
        write_csv(release / f'{role}_cases.csv', cases)
        files[f'{role}_cases'] = f'{role}_cases.csv'
        baseline = rng.uniform(.1, 1., (count, N_GENES)).astype('float32')
        np.save(release / f'{role}_baseline.npy', baseline)
        files[f'{role}_baseline'] = f'{role}_baseline.npy'
        if role != 'test':
            if delta_value is not None:
                delta = np.full((count, N_GENES), delta_value[role], dtype='float32')
            else:
                delta = rng.uniform(-.05, .05, (count, N_GENES)).astype('float32')
            np.save(release / f'{role}_delta.npy', delta)
            files[f'{role}_delta'] = f'{role}_delta.npy'
    manifest = {
        'schema': 'stage2-02-models-tenv-inputs-v1',
        'target_kind': 'delta',
        'n_genes': N_GENES,
        'declared': {'folds': 2, 'seeds': 2},
        'identity': {'fold_row': 0, 'fold': 0, 'seed_row': 0, 'seed': 17},
        'roles': {role: count for role, _, count in roles},
        'morgan': {'radius': 2, 'bits': 2048, 'chirality': True, 'binary': True},
        'files': files,
        'sha256': {name: sha256(release / relative) for name, relative in files.items()},
        'source': {'corpus_sha256': 'synthetic', 'folds_sha256': 'synthetic',
                   'seeds_sha256': 'synthetic', 'cohorts_sha256': 'synthetic'},
    }
    (release / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    request = {
        'identity': {'model': 'prnet', **manifest['identity']},
        'release': str(release),
        'n_genes': N_GENES,
        'target_kind': 'delta',
        'output_directory': str(root / 'output'),
    }
    (root / 'request.json').write_text(json.dumps(request, indent=2) + '\n')
    return request


def test_environment(temporary):
    environment = dict(
        OMP_NUM_THREADS='2', MKL_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2',
        NUMEXPR_NUM_THREADS='2',
        NUMBA_CACHE_DIR=str(Path(temporary) / 'numba'),
        MPLCONFIGDIR=str(Path(temporary) / 'matplotlib'),
        PATH='/usr/local/bin:/usr/bin:/bin',
    )
    return environment


def run_agent(root, config_path, extra=()):
    return subprocess.run(
        [PYTHON, str(ADAPTER), '--request', str(root / 'request.json'),
         '--config', str(config_path), *extra],
        capture_output=True, text=True, cwd=str(REPO),
        env=test_environment(root))


def bounded_config(root, *, max_epochs=2, max_training_seconds=600,
                   no_strict_decrease_epochs=5, learning_rate=None,
                   inference_max_seconds=300):
    config = json.loads(BASE_CONFIG.read_text())
    config.update(max_epochs=max_epochs, batch_size=2, device='cpu',
                  no_strict_decrease_epochs=no_strict_decrease_epochs,
                  above_running_min_epochs=20,
                  max_training_seconds=max_training_seconds,
                  inference_max_seconds=inference_max_seconds)
    if learning_rate is not None:
        config['learning_rate'] = learning_rate
    config['early_stopping'].update(
        max_epochs=max_epochs,
        no_strict_decrease_epochs=config['no_strict_decrease_epochs'],
        above_running_min_epochs=config['above_running_min_epochs'])
    path = root / 'config.json'
    path.write_text(json.dumps(config, indent=2) + '\n')
    return path


class CompletedRun(unittest.TestCase):
    def test_trained_checkpoint_reloads_and_predicts_the_declared_shape(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rng = np.random.default_rng(44)
            request = write_release(root, rng)
            config_path = bounded_config(root)
            process = run_agent(root, config_path)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            output = Path(request['output_directory'])

            prediction = np.load(output / 'prediction.npy')
            self.assertEqual(prediction.shape, (2, N_GENES))
            self.assertEqual(prediction.dtype, np.float32)
            self.assertTrue(np.isfinite(prediction).all())

            self.assertEqual(json.loads((output / 'case_ids.json').read_text()),
                             ['test0', 'test1'])
            checkpoint = json.loads((output / 'checkpoint.json').read_text())
            self.assertEqual(checkpoint['identity'], request['identity'])
            self.assertEqual(checkpoint['author_objective'], 'GaussianNLLLoss')
            self.assertEqual(checkpoint['adapter_parameters'], 2048 * 1024 + 1024)
            self.assertEqual(checkpoint['n_genes'], N_GENES)
            self.assertEqual(checkpoint['n_cases'], 2)
            self.assertEqual(checkpoint['selection'],
                             'minimum author Gaussian NLL on validation; '
                             'fixed validation noise stream')
            self.assertEqual(checkpoint['behavior']['reload_max_abs_error'], 0.)
            for key in ('adapter_parameter_change_l2',
                        'predictor_parameter_change_l2',
                        'baseline_sensitivity_l2'):
                self.assertGreater(checkpoint['behavior'][key], 0, key)
            self.assertTrue(np.isfinite(checkpoint['best_validation_loss']))
            self.assertEqual(len(checkpoint['history']), 2)
            self.assertEqual(checkpoint['stop_reason'], 'max_epochs')
            self.assertEqual(checkpoint['config']['max_epochs'], 2)
            self.assertGreater(checkpoint['peak_rss_kbytes'], 0)

            stopping = checkpoint['early_stopping']
            self.assertEqual(stopping['rules']['no_strict_decrease_epochs'], 5)
            self.assertEqual(stopping['rules']['above_running_min_epochs'], 20)
            self.assertEqual(stopping['rules']['max_epochs'], 2)
            self.assertIn('5-epoch rule dominates', stopping['dominance'])
            self.assertTrue(stopping['dominance_observed'])

            terminal = json.loads((output / 'terminal.json').read_text())
            self.assertEqual(terminal['status'], 'complete')
            self.assertEqual(terminal['phase'], 'complete')
            self.assertIsNone(terminal['error'])
            self.assertFalse((output / 'prediction.npy.partial').exists())

            identity = json.loads((output / 'identity.json').read_text())
            self.assertEqual(identity['identity'], request['identity'])
            self.assertEqual(identity['n_genes'], N_GENES)
            self.assertEqual(identity['n_cases'], 2)
            self.assertEqual(identity['prediction_sha256'],
                             sha256(output / 'prediction.npy'))
            self.assertEqual(identity['checkpoint_sha256'],
                             sha256(output / 'best.pt'))
            self.assertEqual(identity['config_sha256'],
                             sha256(output / 'config.json'))
            self.assertEqual(identity['release_manifest_sha256'],
                             sha256(Path(request['release']) / 'manifest.json'))
            self.assertTrue(identity['code_sha256'])
            self.assertEqual(identity['environment']['executable'], PYTHON)
            self.assertTrue(identity['environment']['upstream_source'].endswith(
                'pert_model/PRnet/source'))
            self.assertEqual(
                len(identity['environment']['upstream_prnet_module_sha256']), 64)

            progress = json.loads((output / 'progress.json').read_text())
            self.assertEqual(progress['phase'], 'complete')
            self.assertGreater(progress['training_seconds'], 0.)

            status = subprocess.run(
                ['git', 'status', '--porcelain', 'pert_model/PRnet/source'],
                capture_output=True, text=True, cwd=str(REPO))
            tracked_changes = [line for line in status.stdout.splitlines()
                               if not line.startswith('??')]
            self.assertEqual(tracked_changes, [],
                             'upstream PRnet source must stay unmodified')


class HonestTerminalStatus(unittest.TestCase):
    def test_training_budget_exhaustion_records_over_budget_without_prediction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rng = np.random.default_rng(7)
            request = write_release(root, rng)
            config_path = bounded_config(root, max_training_seconds=0)
            process = run_agent(root, config_path)
            self.assertNotEqual(process.returncode, 0)
            output = Path(request['output_directory'])
            self.assertFalse((output / 'prediction.npy').exists())
            self.assertFalse((output / 'prediction.npy.partial').exists())
            terminal = json.loads((output / 'terminal.json').read_text())
            self.assertEqual(terminal['status'], 'over_budget')
            self.assertEqual(terminal['phase'], 'training')
            self.assertIn('budget', terminal['error'])
            self.assertIsNone(terminal.get('prediction'))

    def test_a_rejected_release_records_failed_instead_of_writing_nothing_silent(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rng = np.random.default_rng(9)
            request = write_release(root, rng)
            manifest_path = Path(request['release']) / 'manifest.json'
            manifest = json.loads(manifest_path.read_text())
            manifest['n_genes'] = 3
            manifest_path.write_text(json.dumps(manifest))
            config_path = bounded_config(root)
            process = run_agent(root, config_path)
            self.assertNotEqual(process.returncode, 0)
            output = Path(request['output_directory'])
            self.assertFalse((output / 'prediction.npy').exists())
            terminal = json.loads((output / 'terminal.json').read_text())
            self.assertEqual(terminal['status'], 'failed')
            self.assertEqual(terminal['phase'], 'validate')
            self.assertIn('n_genes', terminal['error'])


class ReloadAndReplay(unittest.TestCase):
    def test_a_best_epoch_before_the_last_is_reloaded_and_replays_exactly(self):
        # Train targets positive, validation targets negative: training improves
        # the fit while worsening validation, so the selected checkpoint is an
        # early epoch and the final in-memory state is not the selected one.
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rng = np.random.default_rng(123)
            request = write_release(
                root, rng, delta_value={'train': .2, 'validation': -.2})
            config_path = bounded_config(
                root, max_epochs=6, no_strict_decrease_epochs=2,
                learning_rate=0.01)
            process = run_agent(root, config_path)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            output = Path(request['output_directory'])
            checkpoint = json.loads((output / 'checkpoint.json').read_text())
            last_epoch = checkpoint['history'][-1]['epoch']
            self.assertEqual(checkpoint['stop_reason'], 'no_strict_decrease')
            self.assertLess(checkpoint['best_epoch'], last_epoch,
                            'this fixture must select an epoch before the last')
            terminal_before = (output / 'terminal.json').read_bytes()
            replay = run_agent(root, config_path, extra=['--replay'])
            self.assertEqual(replay.returncode, 0, replay.stdout + replay.stderr)
            published = np.load(output / 'prediction.npy')
            replayed = np.load(output / 'replay.npy')
            np.testing.assert_array_equal(replayed, published)
            summary = json.loads(replay.stdout.strip().splitlines()[-1])
            self.assertTrue(summary['matches_published'])
            self.assertEqual(summary['replay_sha256'], sha256(output / 'replay.npy'))
            self.assertEqual((output / 'terminal.json').read_bytes(), terminal_before,
                             'replay must never rewrite the run terminal record')
            self.assertLessEqual(checkpoint['best_epoch'], 6)


class ConfigContract(unittest.TestCase):
    def test_epoch_ceiling_above_the_adr_limit_is_rejected_before_training(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            request = write_release(root, np.random.default_rng(5))
            config_path = bounded_config(root, max_epochs=101)
            process = run_agent(root, config_path)
            self.assertNotEqual(process.returncode, 0)
            output = Path(request['output_directory'])
            self.assertFalse((output / 'prediction.npy').exists())
            terminal = json.loads((output / 'terminal.json').read_text())
            self.assertEqual(terminal['status'], 'failed')
            self.assertEqual(terminal['phase'], 'validate')
            self.assertIn('max_epochs', terminal['error'])

    def test_inference_budget_exhaustion_records_over_budget_after_training(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            request = write_release(root, np.random.default_rng(11))
            config_path = bounded_config(root, max_epochs=2,
                                         inference_max_seconds=1e-9)
            process = run_agent(root, config_path)
            self.assertNotEqual(process.returncode, 0)
            output = Path(request['output_directory'])
            self.assertFalse((output / 'prediction.npy').exists())
            self.assertFalse((output / 'prediction.npy.partial').exists())
            terminal = json.loads((output / 'terminal.json').read_text())
            self.assertEqual(terminal['status'], 'over_budget')
            self.assertEqual(terminal['phase'], 'inference')
            self.assertIn('inference budget', terminal['error'])
            # Training completed, so the terminal record carries the observed
            # early-stopping accounting (the ADR 0063 dominance evidence).
            self.assertIsNotNone(terminal['early_stopping'])
            self.assertEqual(
                terminal['early_stopping']['rules']['no_strict_decrease_epochs'], 5)
            self.assertIn('5-epoch rule dominates',
                          terminal['early_stopping']['dominance'])


if __name__ == '__main__':
    unittest.main()
