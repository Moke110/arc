"""Replay the frozen Step 04 evaluator CLI with small artifacts for both datasets."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


CLI = Path(__file__).resolve().parents[1] / 'script/evaluate.py'
METRICS = ('AUPRC_positive@100', 'AUPRC_negative@100', 'Recall_positive@P80',
           'Recall_negative@P80', 'PCC', 'Spearman', 'RMSE', 'R2')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + '\n')
    return path


def csv_rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def release(root, dataset):
    fold = root / 'fold_0'
    evaluator = fold / 'evaluator'
    evaluator.mkdir(parents=True)
    (root / 'genes.csv').write_text('gene_id\na\nb\nc\nd\n')
    (root / 'folds.csv').write_text('frozen split\n')
    header = 'case_id,group_id,drug_key,dose_um,control_index\n'
    (fold / 'train.csv').write_text(header + 't1,t,train,1,0\nt2,t,train,1,0\nu1,u,train2,1,0\n')
    # Training group means are [2, 0, -2, 0] and [0, 2, 0, -2].
    np.save(fold / 'train_expression.npy', np.array([[4, 0, -4, 0], [0, 0, 0, 0],
                                                    [0, 2, 0, -2]], dtype='float32'))
    (fold / 'test_cases.csv').write_text(header + 'x1,x,test,1,0\ny1,y,test,2,1\nx2,x,test,1,2\n')
    (evaluator / 'groups.csv').write_text('group_id,n_wells,drug_key,dose_um,cell_line_id\nx,2,test,1,L\ny,1,test,2,L\n')
    # Step 03's perfect and undefined sanity fixture.
    np.save(evaluator / 'truth.npy', np.array([[2, 1, -1, -2], [0, 0, 0, 0]], dtype='float32'))
    np.save(evaluator / 'positive.npy', [[True, True, False, False], [False] * 4])
    np.save(evaluator / 'negative.npy', [[False, False, True, True], [False] * 4])
    paths = {name: str(path.relative_to(root)) for name, path in {
        'genes': root / 'genes.csv', 'train_cases': fold / 'train.csv',
        'train_expression': fold / 'train_expression.npy', 'test_cases': fold / 'test_cases.csv',
    }.items()}
    files = {str(path.relative_to(root)): digest(path) for path in root.rglob('*') if path.is_file()}
    manifest = write(root / 'manifest.json', {'schema': f'{dataset}-input-v1', 'target_kind': 'delta',
                                              'gene_count': 4, 'files': files})
    consumer = {'schema': 'step04-input-v1', 'dataset': dataset, 'fold': 0, 'n_genes': 4,
                'target_kind': 'delta', 'source_input_sha256': digest(manifest), 'corpus_sha256': {},
                'preprocessing': {}, 'split_sha256': files['folds.csv'], 'gene_sha256': files['genes.csv'],
                'predictor_files': paths, 'file_sha256': {name: files[path] for name, path in paths.items()}}
    identity = {key: consumer[key] for key in ('dataset', 'fold', 'source_input_sha256', 'split_sha256', 'gene_sha256')}
    identity['input_sha256'] = json_hash(consumer)
    request = write(root / 'request.json', {'schema': 'step04-predictor-v1', 'identity': identity,
                                            'input_root': str(root), 'n_genes': 4, 'target_kind': 'delta',
                                            'predictor_files': {k: str(root / p) for k, p in paths.items()},
                                            'file_sha256': consumer['file_sha256'], 'input': consumer})
    return request, identity


def prediction(directory, identity, seed):
    directory.mkdir()
    identity = {**identity, 'seed': seed, 'stage': 'fixture'}
    recipe = {'stage': 'fixture', 'seed': seed, 'fold': identity['fold']}
    recipe_sha = json_hash(recipe)
    # Freeze opaque bytes representing the trainer's persisted checkpoint.
    checkpoint = directory / 'best.pt'
    checkpoint.write_bytes(f'fixture checkpoint {identity}'.encode())
    contract = {'identity': identity, 'recipe_sha256': recipe_sha, 'best_step': 2,
                'checkpoint_sha256': digest(checkpoint), 'best_validation_loss': .5,
                'model_state_sha256': hashlib.sha256(b'fixture weights').hexdigest()}
    run = write(directory / 'run.json', {**contract, 'recipe': recipe, 'status': 'completed',
                                         'monitor': 'validation_validation_loss_epoch',
                                         'history': [{'step': 1, 'validation_loss': 1.0},
                                                     {'step': 2, 'validation_loss': .5}]})
    np.save(directory / 'prediction.npy', np.array([[4, 2, -2, -4], [1, 1, 1, 1], [0, 0, 0, 0]],
                                                    dtype='float32'))
    cases = write(directory / 'cases.json', [['x1', 'x'], ['y1', 'y'], ['x2', 'x']])
    manifest = write(directory / 'prediction.json', {
        'schema': 'step04-heldout-prediction-v1', **contract,
        'run_sha256': digest(run), 'cases_sha256': digest(cases),
        'prediction_sha256': digest(directory / 'prediction.npy'),
    })
    return recipe_sha, manifest


class EvaluatorReplay(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.folds = []
        self.outcomes = []
        for dataset in ('tahoe', 'tenv'):
            request, identity = release(self.root / dataset, dataset)
            recipe_sha, manifest = prediction(self.root / f'{dataset}_run', identity, 17)
            self.folds.append({'request': str(request.relative_to(self.root)), 'request_sha256': digest(request),
                               'recipe_sha256_by_seed': {'17': recipe_sha}})
            self.outcomes.append({'dataset': dataset, 'fold': 0, 'seed': 17, 'status': 'completed',
                                  'prediction': str(manifest.relative_to(self.root)),
                                  'prediction_manifest_sha256': digest(manifest)})

    def evaluate(self, destination):
        plan = write(self.root / 'plan.json', {'schema': 'step04-evaluation-plan-v1',
                                               'stage': 'fixture', 'folds': self.folds})
        outcomes = write(self.root / 'outcomes.json', {'schema': 'step04-evaluation-outcomes-v1',
                                                       'outcomes': self.outcomes})
        return subprocess.run([sys.executable, str(CLI), str(plan), str(outcomes), str(destination)],
                              capture_output=True, text=True, timeout=30)

    def test_matches_references_and_replays_all_eight_metrics(self):
        output = self.root / 'scores'
        result = self.evaluate(output)
        self.assertEqual(result.returncode, 0, result.stderr)
        runs = csv_rows(output / 'per_run.csv')
        for dataset in ('tahoe', 'tenv'):
            references = output / 'references' / f'{dataset}_fold_0'
            self.assertEqual(json.loads((references / 'group_ids.json').read_text()), ['x', 'y'])
            np.testing.assert_array_equal(np.load(references / 'no_change.npy'), [0, 0, 0, 0])
            np.testing.assert_array_equal(np.load(references / 'delta_mean.npy'), [1, 1, -1, -1])
            dataset_runs = {row['model']: row for row in runs if row['dataset'] == dataset}
            self.assertEqual(set(dataset_runs), {'perturbdiff', 'no_change', 'delta_mean'})
            model = dataset_runs['perturbdiff']
            for metric in METRICS:
                self.assertAlmostEqual(float(model[metric]), .5 if metric == 'RMSE' else 1.)
                self.assertEqual(int(model[f'{metric}_valid_groups']), 2 if metric == 'RMSE' else 1)
            no_change = dataset_runs['no_change']
            self.assertEqual(no_change['PCC'], '')
            self.assertEqual(no_change['PCC_valid_groups'], '0')
            self.assertEqual(no_change['AUPRC_positive@100'], '0.5')
        self.assertEqual(len(runs), 6)
        replay = self.root / 'replay'
        result = self.evaluate(replay)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual({path.relative_to(output): digest(path) for path in output.rglob('*') if path.is_file()},
                         {path.relative_to(replay): digest(path) for path in replay.rglob('*') if path.is_file()})


    def test_preserves_missing_and_failed_slots_without_reference_seed_copies(self):
        self.folds[0]['recipe_sha256_by_seed'].update({'42': 'a' * 64, '2026': 'b' * 64})
        self.outcomes.append({'dataset': 'tahoe', 'fold': 0, 'seed': 42,
                              'status': 'failed', 'reason': 'worker_error'})
        output = self.root / 'partial'
        result = self.evaluate(output)
        self.assertEqual(result.returncode, 0, result.stderr)
        tahoe = [row for row in csv_rows(output / 'per_run.csv') if row['dataset'] == 'tahoe']
        self.assertEqual([(row['model'], row['seed'], row['status']) for row in tahoe], [
            ('perturbdiff', '17', 'completed'), ('perturbdiff', '42', 'failed'),
            ('perturbdiff', '2026', 'missing'), ('no_change', '', 'completed'),
            ('delta_mean', '', 'completed'),
        ])
        self.assertEqual(tahoe[1]['reason'], 'worker_error')
        self.assertEqual(tahoe[2]['RMSE'], '')
        self.assertEqual(tahoe[2]['RMSE_valid_groups'], '0')
        summary = csv_rows(output / 'summary.csv')[0]
        self.assertEqual((summary['expected_runs'], summary['completed_runs'], summary['RMSE_valid_runs']),
                         ('3', '1', '1'))
        audit = json.loads((output / 'evaluation.json').read_text())
        self.assertFalse(audit['complete'])

    def test_rejects_duplicate_slot(self):
        self.outcomes.append(self.outcomes[0].copy())
        output = self.root / 'invalid'
        result = self.evaluate(output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('duplicate outcome slot', result.stderr)
        self.assertFalse(output.exists())

    def test_rejects_rehashed_prediction_for_different_checkpoint_or_group(self):
        outcome = self.outcomes[0]
        path = self.root / outcome['prediction']
        original = json.loads(path.read_text())
        for failure in ('checkpoint', 'selection', 'group'):
            with self.subTest(failure=failure):
                manifest = original.copy()
                if failure == 'checkpoint':
                    manifest['checkpoint_sha256'] = 'a' * 64
                    message = 'checkpoint identity'
                elif failure == 'selection':
                    run_path = path.parent / 'run.json'
                    run_original = run_path.read_text()
                    run = json.loads(run_original)
                    run.update(best_step=1, best_validation_loss=1.)
                    manifest.update(best_step=1, best_validation_loss=1.,
                                    run_sha256=digest(write(run_path, run)))
                    message = 'strict validation minimum'
                else:
                    cases = write(path.parent / 'cases.json', [['x1', 'y'], ['y1', 'x'], ['x2', 'x']])
                    manifest['cases_sha256'] = digest(cases)
                    message = 'case and group identities'
                outcome['prediction_manifest_sha256'] = digest(write(path, manifest))
                output = self.root / 'invalid'
                result = self.evaluate(output)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
                self.assertFalse(output.exists())
                if failure == 'selection':
                    run_path.write_text(run_original)

    def test_rejects_frozen_evaluator_group_count_mismatch(self):
        fold = self.folds[0]
        request_path = self.root / fold['request']
        request = json.loads(request_path.read_text())
        source = Path(request['input_root'])
        groups_path = source / 'fold_0/evaluator/groups.csv'
        groups_path.write_text(groups_path.read_text().replace('x,2,', 'x,1,'))
        manifest_path = source / 'manifest.json'
        manifest = json.loads(manifest_path.read_text())
        manifest['files']['fold_0/evaluator/groups.csv'] = digest(groups_path)
        source_hash = digest(write(manifest_path, manifest))
        request['input']['source_input_sha256'] = source_hash
        request['identity'].update(source_input_sha256=source_hash, input_sha256=json_hash(request['input']))
        fold['request_sha256'] = digest(write(request_path, request))
        self.outcomes.pop(0)
        output = self.root / 'invalid'
        result = self.evaluate(output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('group membership differs from test cases', result.stderr)
        self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
