"""Exercise warm-start split and prediction gates at the artifact CLI."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from collections import Counter

import numpy as np

STAGE2 = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists()) / 'stage2'
STEP05 = STAGE2 / '05_pdiff_transfer'
STEP03 = STAGE2 / '03_modes_tenv'
STEP04 = STAGE2 / '04_pdiff_hyperparameter'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def read_rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_rows(path, cases):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)


class WarmTenvCLI(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.formal = self.root / 'formal'
        self.source = self.root / 'source'
        self.output = self.root / 'warm'
        self.score = self.root / 'score'
        self.runtime = self.root / 'runtime'
        write(self.runtime / 'stage2/03_modes_tenv/execution/input_record/verification.json',
              {'manifest_sha256': self.make_release()})
        request = self.formal / 'request.json'
        write(request, {'input_root': str(self.source)})
        recipes = {}
        for seed in (17, 42, 2026):
            recipe = self.formal / f'seed-{seed}.json'
            write(recipe, {'fold': 0, 'stage': 'formal', 'seed': seed, 'max_steps': 30000})
            recipes[str(seed)] = {'path': recipe.name, 'sha256': sha(recipe)}
        plan = {'schema': 'step04-evaluation-plan-v1', 'stage': 'formal', 'dataset': 'tenv',
                'folds': [{'request': request.name, 'request_sha256': sha(request), 'recipes': recipes,
                           'recipe_sha256_by_seed': dict.fromkeys(recipes, 'a' * 64)}],
                'inference': {}}
        for root in (self.formal, self.runtime / 'stage2/04_pdiff_hyperparameter/execution/formal_record/tenv'):
            write(root / 'plan.json', plan)

    def make_release(self):
        source = self.source
        source.mkdir()
        drugs = [{'drug_key': f'D{n}', 'parent_connectivity': f'P{n // 2}', 'morgan_index': n}
                 for n in range(24)]
        self.wells = [
            {'case_id': f'{drug}-{dose}-{repeat}', 'drug_key': drug, 'dose_um': dose,
             'control_index': '0', 'group_id': f'{drug}-{dose}'}
            for drug in (row['drug_key'] for row in drugs) for dose in ('0.2', '1.0', '10.0')
            for repeat in (1, 2)]
        write_rows(source / 'drugs.csv', drugs)
        write_rows(source / 'wells.csv', self.wells)
        (source / 'genes.csv').write_text('symbol\n' + ''.join(f'gene{n}\n' for n in range(18164)))
        np.save(source / 'controls.npy', np.zeros((1, 18164), dtype='float32'))
        np.save(source / 'morgan.npy', np.zeros((24, 2048), dtype='float32'))
        for fold in range(3):
            for role, wanted in (('train', range(0, 96)), ('validation', range(96, 144))):
                selected = [self.wells[i] for i in wanted]
                path = source / f'fold_{fold}' / f'{role}.csv'
                write_rows(path, selected)
                values = np.zeros((len(selected), 18164), dtype='float32')
                values[:, 0] = [index + 1 for index in wanted]
                values[:, 1] = -values[:, 0]
                np.save(path.with_name(f'{role}_expression.npy'), values)
        files = {str(path.relative_to(source)): sha(path) for path in source.rglob('*') if path.is_file()}
        manifest = source / 'manifest.json'
        write(manifest, {'schema': 'tenv-input-v1', 'target_kind': 'delta', 'gene_count': 18164, 'files': files})
        return sha(manifest)

    def cli(self, *args, succeeds=True):
        result = subprocess.run([sys.executable, str(STEP05 / 'execution/script/warm_tenv.py'),
                                 '--runtime', str(self.runtime), *map(str, args)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, succeeds, result.stderr)
        return result.stderr

    def test_prepares_balanced_atomic_split_and_rejects_resealed_case_leakage(self):
        self.cli('prepare', self.formal, self.output, '--smoke')
        plan = self.output / 'plan.json'
        self.cli('verify', self.output, '--plan-sha256', sha(plan))
        split = json.loads((self.output / 'split.json').read_text())
        roles = {(item['parent_connectivity'], item['dose_um']): item['role'] for item in split['combinations']}
        self.assertEqual(len(roles), 36)
        self.assertEqual(Counter(dose for (_, dose), role in roles.items() if role == 'test'),
                         {'0.2': 4, '1': 4, '10': 4})
        self.assertEqual(sum(role == 'validation' for role in roles.values()), 4)
        self.assertEqual({role for role in roles.values()}, {'train', 'validation', 'test'})
        request_path = self.output / 'request.json'
        request = json.loads(request_path.read_text())
        self.assertNotIn('test_expression', request['predictor_files'])
        self.assertEqual(set(json.loads(plan.read_text())['formal_recipes']), {'17', '42', '2026'})
        test_path = Path(request['predictor_files']['test_cases'])
        cases = read_rows(test_path)
        cases.append(read_rows(request['predictor_files']['train_cases'])[0])
        write_rows(test_path, cases)
        request['file_sha256']['test_cases'] = sha(test_path)
        request['input']['file_sha256']['test_cases'] = sha(test_path)
        request['identity']['input_sha256'] = json_hash(request['input'])
        write(request_path, request)
        plan_doc = json.loads(plan.read_text())
        plan_doc['request_sha256'] = sha(request_path)
        write(plan, plan_doc)
        self.assertIn('test cases differ', self.cli('verify', self.output, '--plan-sha256', sha(plan), succeeds=False))

    def test_score_rejects_absent_prediction_before_opening_source_expression(self):
        self.cli('prepare', self.formal, self.output, '--smoke')
        (self.source / 'fold_0/train_expression.npy').unlink()
        error = self.cli('score', self.output, self.score, '--plan-sha256',
                         sha(self.output / 'plan.json'), succeeds=False)
        self.assertIn('prediction.json', error)
        self.assertFalse(self.score.exists())

    def test_scores_a_sealed_prediction_and_replays_both_references(self):
        self.cli('prepare', self.formal, self.output, '--smoke')
        plan_sha = sha(self.output / 'plan.json')
        request = json.loads((self.output / 'request.json').read_text())
        recipe = json.loads((self.output / 'smoke.json').read_text())
        cases = read_rows(request['predictor_files']['test_cases'])
        fit = self.output / 'slots/seed-17/fit'
        pred = fit.parent / 'prediction'
        pred.mkdir(parents=True)
        (pred / 'best.pt').write_bytes(b'fixture EMA weights')
        contract = {'identity': {**request['identity'], 'stage': 'smoke', 'seed': 17},
                    'recipe_sha256': json_hash(recipe), 'best_step': 6,
                    'best_validation_loss': .5, 'checkpoint_sha256': sha(pred / 'best.pt'),
                    'model_state_sha256': hashlib.sha256(b'weights').hexdigest()}
        reload = {'validation_loss': .5, 'prediction_max_abs_error': 0}
        run = {**contract, 'recipe': recipe, 'status': 'completed', 'reload': reload,
               'monitor': 'validation_validation_loss_epoch',
               'history': [{'step': 6, 'validation_loss': .5}, {'step': 12, 'validation_loss': .75}]}
        write(fit / 'run.json', run)
        write(pred / 'run.json', run)
        write(fit / 'fit.json', {'plan_sha256': plan_sha, 'run_sha256': sha(fit / 'run.json')})
        write(pred / 'replay.json', {**contract, **reload, 'plan_sha256': plan_sha})
        write(pred / 'cases.json', [[case['case_id'], case['group_id']] for case in cases])
        indexes = {row['case_id']: index + 1 for index, row in enumerate(self.wells)}
        exact = np.zeros((len(cases), 18164), dtype='float32')
        exact[:, 0] = [indexes[case['case_id']] for case in cases]
        exact[:, 1] = -exact[:, 0]
        np.save(pred / 'prediction.npy', exact)
        write(pred / 'prediction.json', {
            **contract, 'schema': 'step04-heldout-prediction-v1', 'plan_sha256': plan_sha,
            'case_count': len(cases), 'group_count': len({case['group_id'] for case in cases}),
            **{f'{key}_sha256': sha(pred / name) for key, name in (
                ('prediction', 'prediction.npy'), ('cases', 'cases.json'),
                ('run', 'run.json'), ('replay', 'replay.json'))},
        })
        self.cli('score', self.output, self.score, '--plan-sha256', plan_sha)
        scores = {row['model']: row for row in read_rows(self.score / 'per_run.csv')}
        self.assertEqual(set(scores), {'perturbdiff', 'no_change', 'delta_mean'})
        self.assertEqual(float(scores['perturbdiff']['RMSE']), 0)
        self.assertEqual(float(scores['perturbdiff']['AUPRC_positive@100']), 1)
        self.assertEqual(float(scores['perturbdiff']['AUPRC_negative@100']), 1)
        self.assertEqual(scores['no_change']['PCC'], '')
        self.assertGreater(float(scores['delta_mean']['RMSE']), 0)
        replay = self.root / 'replay'
        self.cli('score', self.output, replay, '--plan-sha256', plan_sha)
        self.assertEqual((self.score / 'evaluation.json').read_bytes(), (replay / 'evaluation.json').read_bytes())


if __name__ == '__main__':
    unittest.main()
