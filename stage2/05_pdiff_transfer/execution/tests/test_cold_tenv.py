"""Exercise Step 05 cold artifact reuse at its CLI boundary."""
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

STEP05 = Path(__file__).resolve().parents[2]
STAGE2 = STEP05.parent
STEP04 = STAGE2 / '04_pdiff_hyperparameter'
STEP03 = STAGE2 / '03_modes_tenv'
sys.path.insert(0, str(STEP04 / 'execution/tests'))
sys.path.insert(0, str(STEP04 / 'execution/script'))
sys.path.insert(0, str(STEP03 / 'execution/script'))

from test_evaluate import csv_rows, digest, json_hash, prediction, write
from evaluate import write_results
from training_inputs import predictor_paths, request_for

SEEDS = (17, 42, 2026)
N_GENES = 18164


class ColdTenvReplay(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        (self.root / '.git').touch()
        self.step04 = self.root / 'stage2' / STEP04.name
        self.record = self.step04 / 'execution/formal_record/tenv'
        self.formal = self.root / 'formal'
        self.checkpoints = self.root / 'checkpoints'
        self.output = self.root / 'output'
        self.cli = self.root / 'stage2' / STEP05.name / 'execution/script/cold_tenv.py'
        self.cli.parent.mkdir(parents=True)
        cli = STEP05 / 'execution/script/cold_tenv.py'
        if cli.is_file():
            shutil.copyfile(cli, self.cli)
        scripts = self.step04 / 'execution/script'
        scripts.mkdir(parents=True)
        for name in ('evaluate.py', 'training_inputs.py', 'training_policy.py'):
            shutil.copyfile(STEP04 / 'execution/script' / name, scripts / name)
        step03 = self.root / 'stage2' / STEP03.name
        (step03 / 'execution').mkdir(parents=True)
        (step03 / 'execution/script').symlink_to(STEP03 / 'execution/script', target_is_directory=True)

        source = self.root / 'input'
        fold = source / 'fold_0'
        evaluator = fold / 'evaluator'
        evaluator.mkdir(parents=True)
        (source / 'genes.csv').write_text('symbol\n' + ''.join(f'gene{n}\n' for n in range(N_GENES)))
        (source / 'drugs.csv').write_text(
            'drug_key,parent_connectivity\ntrain,P1\ntrain2,P2\nvalidation,P3\ntest,P4\n')
        (source / 'folds.csv').write_text('fold,drug_key,parent_connectivity,role\n'
                                        '0,train,P1,train\n0,train2,P2,train\n'
                                        '0,validation,P3,validation\n0,test,P4,test\n')
        (source / 'controls.csv').write_text('control_index,plate,n_wells\n0,L,2\n1,R,2\n2,S,2\n')
        np.save(source / 'controls.npy', np.zeros((3, N_GENES), dtype='float32'))
        np.save(source / 'morgan.npy', np.zeros((4, 2048), dtype='float32'))
        header = 'case_id,group_id,drug_key,dose_um,control_index\n'
        train = 't1,t,train,1,0\nt2,t,train,1,0\nu1,u,train2,1,0\n'
        validation = 'v1,v,validation,1,0\n'
        test = 'x1,x,test,1,0\ny1,y,test,2,1\nx2,x,test,1,2\n'
        (source / 'wells.csv').write_text(header + train + validation + test)
        for filename, memberships in (('train', train), ('validation', validation), ('test_cases', test)):
            (fold / f'{filename}.csv').write_text(header + memberships)
        np.save(fold / 'train_expression.npy', full_axis([[4, 0, -4, 0], [0, 0, 0, 0], [0, 2, 0, -2]]))
        np.save(fold / 'validation_expression.npy', full_axis([[0, 1, 0, -1]]))
        (evaluator / 'groups.csv').write_text('group_id,n_wells,drug_key,dose_um\nx,2,test,1\ny,1,test,2\n')
        np.save(evaluator / 'truth.npy', full_axis([[2, 1, -1, -2], [0, 0, 0, 0]]))
        np.save(evaluator / 'positive.npy', full_axis([[1, 1, 0, 0], [0, 0, 0, 0]]).astype(bool))
        np.save(evaluator / 'negative.npy', full_axis([[0, 0, 1, 1], [0, 0, 0, 0]]).astype(bool))
        paths = predictor_paths(0)
        files = {str(path.relative_to(source)): digest(path) for path in source.rglob('*') if path.is_file()}
        manifest = write(source / 'manifest.json', {'schema': 'tenv-input-v1', 'target_kind': 'delta',
                                                    'gene_count': N_GENES, 'files': files})
        write(step03 / 'execution/input_record/verification.json', {'manifest_sha256': digest(manifest)})
        consumer = {'schema': 'step04-input-v1', 'dataset': 'tenv', 'fold': 0, 'n_genes': N_GENES,
                    'target_kind': 'delta', 'source_input_sha256': digest(manifest), 'corpus_sha256': None,
                    'preprocessing': None, 'split_sha256': files['folds.csv'], 'gene_sha256': files['genes.csv'],
                    'predictor_files': paths, 'file_sha256': {name: files[path] for name, path in paths.items()}}
        request = request_for(source, consumer)
        request_path = write(self.formal / 'requests/fold-0.json', request)
        recipes = {str(seed): {'stage': 'formal', 'seed': seed, 'fold': 0} for seed in SEEDS}
        selection = write(self.formal / 'selection.json', {})
        entry = {'request': 'requests/fold-0.json', 'request_sha256': digest(request_path),
                 'policy_sha256': '1' * 64, 'recipes': {},
                 'recipe_sha256_by_seed': {seed: json_hash(recipe) for seed, recipe in recipes.items()}}
        for seed, recipe in recipes.items():
            name = f'recipes/fold-0/seed-{seed}.json'
            entry['recipes'][seed] = {'path': name, 'sha256': digest(write(self.formal / name, recipe))}
        plan = write(self.formal / 'plan.json', {'schema': 'step04-evaluation-plan-v1', 'stage': 'formal',
                                                'dataset': 'tenv', 'selection_sha256': digest(selection),
                                                'folds': [entry, {}, {}]})
        plan_sha = digest(plan)
        inventory, index, arrays = [], [], {}
        for seed in SEEDS:
            local = self.checkpoints / f'fold-0/seed-{seed}'
            artifact = self.formal / f'slots/fold-0-seed-{seed}/prediction'
            artifact.parent.mkdir(parents=True)
            _, path = prediction(artifact, request['identity'], seed, 'formal', plan_sha)
            metadata = json.loads(path.read_text())
            values = full_axis([[4, 2, -2, -4], [1, 1, 1, 1], [0, 0, 0, 0]])
            np.save(artifact / 'prediction.npy', values)
            metadata.update(prediction_sha256=digest(artifact / 'prediction.npy'), case_count=3, group_count=2)
            write(path, metadata)
            local.mkdir(parents=True)
            for name in ('best.pt', 'run.json', 'replay.json'):
                shutil.copyfile(artifact / name, local / name)
            shutil.copyfile(request_path, local / 'request.json')
            recipe_entry = entry['recipes'][str(seed)]
            shutil.copyfile(self.formal / recipe_entry['path'], local / 'recipe.json')
            index.append({key: metadata[key] for key in ('identity', 'recipe_sha256', 'checkpoint_sha256',
                                                         'model_state_sha256', 'best_step', 'best_validation_loss',
                                                         'run_sha256', 'replay_sha256')} | {
                'request_sha256': digest(request_path), 'policy_sha256': entry['policy_sha256'],
                'raw_recipe_sha256': recipe_entry['sha256'], 'directory': f'fold-0/seed-{seed}'})
            inventory.append({'dataset': 'tenv', 'fold': 0, 'seed': seed, 'status': 'completed',
                              'prediction': str(path.relative_to(self.formal)),
                              'prediction_manifest_sha256': digest(path)})
            arrays[('tenv', 0, seed)] = values
        outcomes = write(self.formal / 'outcomes.json', {'schema': 'step04-evaluation-outcomes-v1',
                                                        'plan_sha256': plan_sha, 'outcomes': inventory})
        self.record.mkdir(parents=True)
        fold_inputs = {'identity': request['identity'], 'root': source, 'files': files, 'paths': paths,
                       'n_genes': N_GENES, 'cases': csv_rows(fold / 'test_cases.csv'),
                       'groups': {'x': [0, 2], 'y': [1]}, 'case_ids': {'x1', 'x2', 'y1'},
                       'recipes': {seed: entry['recipe_sha256_by_seed'][str(seed)] for seed in SEEDS}}
        write_results(self.record, {('tenv', 0): fold_inputs},
                      {('tenv', 0, row['seed']): row for row in inventory}, arrays,
                      {'schema': 'step04-evaluation-v1', 'stage': 'formal', 'datasets': ['tenv'],
                       'plan_sha256': plan_sha, 'outcomes_sha256': digest(outcomes)}, include_pr=True)
        (self.record / 'tenv/pr_curves.json').rename(self.record / 'pr_curves.json')
        write(self.checkpoints / 'checkpoints.json', {'schema': 'step04-tenv-checkpoints-v1',
                                                     'plan_sha256': plan_sha, 'selection_sha256': digest(selection),
                                                     'outcomes_sha256': digest(outcomes),
                                                     'evaluation_sha256': digest(self.record / 'evaluation.json'),
                                                     'checkpoints': index})
        for file in (plan, outcomes, self.checkpoints / 'checkpoints.json'):
            shutil.copyfile(file, self.record / file.name)

    def replay(self):
        return subprocess.run([sys.executable, str(self.cli), str(self.formal), str(self.checkpoints),
                               str(self.output)], capture_output=True, text=True, timeout=30)

    def test_cold_regimes_share_one_verified_cm_and_reference_comparison(self):
        result = self.replay()
        self.assertEqual(result.returncode, 0, result.stderr)
        split = json.loads((self.output / 'split.json').read_text())
        self.assertEqual(split['roles']['test']['parent_connectivity'], ['P4'])
        self.assertEqual(split['roles']['validation']['case_count'], 1)
        comparison = json.loads((self.output / 'comparison.json').read_text())
        self.assertEqual(comparison['regimes']['cold-drug'], comparison['regimes']['cold-cm-drug'])
        self.assertEqual(set(comparison['cm']), {'17', '42', '2026'})
        np.testing.assert_array_equal(np.load(self.output / 'references/tenv_fold_0/delta_mean.npy')[:4],
                                      [1, 1, -1, -1])
        with (self.output / 'per_run.csv').open() as stream:
            result = list(csv.DictReader(stream))
        self.assertEqual([row['seed'] for row in result], ['17', '42', '2026', '', ''])
        self.assertEqual([row['AUPRC_positive@100'] for row in result[:3]], ['1.0'] * 3)

    def test_rejects_checkpoint_from_wrong_seed(self):
        shutil.copyfile(self.checkpoints / 'fold-0/seed-42/best.pt',
                        self.checkpoints / 'fold-0/seed-17/best.pt')
        result = self.replay()
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn('seed-17/best.pt', result.stderr)
        self.assertFalse(self.output.exists())

    def test_rejects_modified_input_manifest(self):
        manifest = self.root / 'input/manifest.json'
        manifest.write_text(manifest.read_text() + '\n')
        result = self.replay()
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn('tenv release differs from the frozen Step 03 identity', result.stderr)
        self.assertFalse(self.output.exists())

    def test_rejects_case_order_even_if_prediction_files_match_their_hashes(self):
        artifact = self.formal / 'slots/fold-0-seed-17/prediction'
        cases = artifact / 'cases.json'
        write(cases, list(reversed(json.loads(cases.read_text()))))
        manifest = artifact / 'prediction.json'
        write(manifest, {**json.loads(manifest.read_text()), 'cases_sha256': digest(cases)})

        outcomes_path = self.formal / 'outcomes.json'
        inventory = json.loads(outcomes_path.read_text())
        inventory['outcomes'][0]['prediction_manifest_sha256'] = digest(manifest)
        write(outcomes_path, inventory)
        shutil.copyfile(outcomes_path, self.record / 'outcomes.json')
        evaluation = self.record / 'evaluation.json'
        write(evaluation, {**json.loads(evaluation.read_text()), 'outcomes_sha256': digest(outcomes_path)})
        index = self.checkpoints / 'checkpoints.json'
        write(index, {**json.loads(index.read_text()), 'outcomes_sha256': digest(outcomes_path),
                      'evaluation_sha256': digest(evaluation)})
        shutil.copyfile(index, self.record / 'checkpoints.json')

        result = self.replay()
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn('prediction case and group identities differ from frozen order', result.stderr)
        self.assertFalse(self.output.exists())


def full_axis(rows):
    return np.pad(np.array(rows, dtype='float32'), ((0, 0), (0, N_GENES - len(rows[0]))))


if __name__ == '__main__':
    unittest.main()
