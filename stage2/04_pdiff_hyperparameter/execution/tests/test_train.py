"""Step 04 artifact CLI monitor, stopping, and input boundary checks."""
import hashlib
import json
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

CLI = Path(__file__).resolve().parents[1] / 'script/train.py'


class TrainingArtifactsTest(unittest.TestCase):
    def test_short_batch_and_strict_best_are_independent_of_step_patience(self):
        with tempfile.TemporaryDirectory() as tmp:
            trace, evidence = (Path(tmp) / name for name in ('trace.json', 'history.json'))
            trace.write_text(json.dumps({
                'minimum_steps': 5, 'patience_steps': 6, 'min_delta': 0.25,
                'events': [
                    {'step': 2, 'batches': [
                        {'loss': [1, 3], 'weights': [0.5, 0.5]},
                        {'loss': [7], 'weights': [1]},
                    ]},
                    {'step': 5, 'batches': [{'loss': [3.875], 'weights': [1]}]},
                    {'step': 8, 'batches': [{'loss': [3.75], 'weights': [1]}]},
                ],
            }))
            subprocess.run([sys.executable, str(CLI), 'replay-monitor', str(trace), str(evidence)],
                           check=True, capture_output=True, text=True)
            report = json.loads(evidence.read_text())
            self.assertEqual(report['monitor'], 'validation_validation_loss_epoch')
            self.assertEqual([event['validation_loss'] for event in report['history']], [4.0, 3.875, 3.75])
            self.assertEqual([event['checkpoint_refreshed'] for event in report['history']], [True] * 3)
            self.assertEqual(report['best_step'], 8)
            self.assertEqual(report['history'][-1]['patience_reference_step'], 2)
            self.assertEqual(report['stop_reason'], 'validation_early_stop')


    def test_frozen_tenv_request_works_without_evaluator_and_rejects_input_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            runtime = Path(tmp)
            release = runtime / 'release'
            fold = release / 'fold_0'
            fold.mkdir(parents=True)
            files = {'genes.csv': b'gene_id\nA\nB\n', 'folds.csv': b'frozen membership',
                     'controls.npy': b'frozen controls', 'morgan.npy': b'frozen fingerprints',
                     'drugs.csv': b'frozen drugs', 'fold_0/train.csv': b'frozen train',
                     'fold_0/validation.csv': b'frozen validation',
                     'fold_0/test_cases.csv': b'frozen test covariates',
                     'fold_0/train_expression.npy': b'frozen train delta',
                     'fold_0/validation_expression.npy': b'frozen validation delta'}
            for name, content in files.items():
                (release / name).write_bytes(content)
            hashes = {name: hashlib.sha256(data).hexdigest() for name, data in files.items()}
            hashes['fold_0/evaluator/truth.npy'] = '0' * 64
            manifest = release / 'manifest.json'
            manifest.write_text(json.dumps({'schema': 'tenv-input-v1', 'target_kind': 'delta',
                                            'gene_count': 18164, 'files': hashes}))
            pin = runtime / 'stage2/03_modes_tenv/execution/input_record/verification.json'
            pin.parent.mkdir(parents=True)
            pin.write_text(json.dumps({'manifest_sha256': hashlib.sha256(manifest.read_bytes()).hexdigest()}))
            request_path = runtime / 'request.json'
            cli = [sys.executable, str(CLI), '--runtime-root', str(runtime)]
            prepared = subprocess.run(cli + ['prepare-tenv', str(release), str(request_path), '--fold', '0'],
                                      capture_output=True, text=True)
            self.assertEqual(prepared.returncode, 0, prepared.stderr)
            self.assertNotIn('evaluator', request_path.read_text())
            request = json.loads(request_path.read_text())
            self.assertEqual(request['identity']['gene_sha256'], hashes['genes.csv'])
            self.assertEqual(request['identity']['split_sha256'], hashes['folds.csv'])
            verified = subprocess.run(cli + ['verify-request', str(request_path)], capture_output=True, text=True)
            self.assertEqual(verified.returncode, 0, verified.stderr)
            vendor = runtime / 'pert_model/PerturbDiff'
            (vendor / 'source').mkdir(parents=True)
            (vendor / 'script').mkdir()
            verifier = CLI.parents[4] / 'pert_model/PerturbDiff/script/verify.sh'
            shutil.copyfile(verifier, vendor / 'script/verify.sh')
            adapter = pin.parent.parent / 'script/perturbdiff_v2_adapter.py'
            adapter.parent.mkdir()
            adapter.write_text("raise RuntimeError('native adapter imported before verification')")
            (runtime / '.git').touch()
            (vendor / 'source/model.py').write_bytes(b'altered')
            pinned = hashlib.sha256(b'original').hexdigest()
            (vendor / 'manifest.sha256').write_text(f'{pinned}  model.py\n')
            fit = subprocess.run(cli + ['fit', str(request_path), 'unused-recipe', 'unused-output'],
                                 capture_output=True, text=True)
            self.assertNotEqual(fit.returncode, 0)
            self.assertIn('vendored PerturbDiff source checksum mismatch', fit.stderr)
            self.assertNotIn('native adapter imported before verification', fit.stderr)
            request['identity']['input_sha256'] = '0' * 64
            bad_request = runtime / 'unbound.json'
            bad_request.write_text(json.dumps(request))
            unbound = subprocess.run(cli + ['verify-request', str(bad_request)], capture_output=True, text=True)
            self.assertNotEqual(unbound.returncode, 0)
            self.assertIn('input identity', unbound.stderr)
            (release / 'genes.csv').write_text('gene_id\nB\nA\n')
            changed = subprocess.run(cli + ['verify-request', str(request_path)], capture_output=True, text=True)
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn('genes', changed.stderr)
            consumer = request['input']
            new_genes = hashlib.sha256((release / 'genes.csv').read_bytes()).hexdigest()
            consumer['gene_sha256'] = request['identity']['gene_sha256'] = new_genes
            consumer['file_sha256']['genes'] = request['file_sha256']['genes'] = new_genes
            request['identity']['input_sha256'] = hashlib.sha256(
                json.dumps(consumer, sort_keys=True, separators=(',', ':')).encode()
            ).hexdigest()
            bad_request.write_text(json.dumps(request))
            forged = subprocess.run(cli + ['verify-request', str(bad_request)], capture_output=True, text=True)
            self.assertNotEqual(forged.returncode, 0)
            self.assertIn('tenv request differs from the frozen Step 03 identity', forged.stderr)


if __name__ == '__main__':
    unittest.main()
