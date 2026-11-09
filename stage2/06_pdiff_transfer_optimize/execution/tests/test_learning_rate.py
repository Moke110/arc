from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

import baseline as base
import learning_rate as lr
from artifacts import write_csv, write_json
from formal_matrix import bind
from training_inputs import json_hash
from training_policy import MONITOR


class LearningRateReplay(unittest.TestCase):
    def test_unavailable_common_step_preserves_native_score_and_seals_incomplete_decision(self):
        with tempfile.TemporaryDirectory() as temp:
            runtime = Path(temp)
            root = runtime / 'learning-rate'
            slot = lr.seed_root(root, 17)
            release = runtime / 'baseline'
            files = {}
            cases = [{'case_id': 'v1', 'group_id': 'g', 'drug_key': 'v',
                      'dose_um': '0.2', 'control_index': '0'}]
            for key in base.INFERENCE_FILES:
                path = runtime / key
                if key == 'validation_cases':
                    write_csv(path, cases)
                else:
                    path.touch()
                files[key] = str(path)
            targets = runtime / 'validation.npy'
            base.np.save(targets, base.np.array([[1, -1]], dtype='float32'))
            files['validation_expression'] = str(targets)
            full = {'identity': {'dataset': 'tenv', 'fold': 0, 'gene_sha256': base.sha256(files['genes'])},
                    'target_kind': 'delta', 'n_genes': 2, 'predictor_files': files,
                    'file_sha256': {key: base.sha256(path) for key, path in files.items()}}
            recipient = runtime / 'request.json'
            write_json(recipient, full)
            request = base.inference_request(full)
            config = {'learning_rate': 0.0002, 'seed': 17, 'stage': 'formal'}
            original, variant = {}, {'learning_rate': 0.0001, 'seed': 17, 'stage': 'formal'}
            for arm, prediction in [('CM', [[1, -1]]), ('Transfer', [[-1, 1]])]:
                path = runtime / f'{arm}.pt'
                path.touch()
                item = {'checkpoint': bind(path), 'metadata': {'best_step': 6500}}
                original[arm] = item
                self.write_prediction(release / arm, item, request, {'arm': arm}, prediction)
            baseline = release / 'plan.json'
            write_json(baseline, {'checkpoints': original, 'request': request, 'case_count': 1})
            baseline_binding = bind(baseline)
            for arm in original:
                path = release / arm / 'prediction.json'
                manifest = base.read_json(path)
                manifest['lineage']['plan_sha256'] = baseline_binding['sha256']
                write_json(path, manifest)

            policy = {'step': 6500, 'cm': 'reuse_selected_ema', 'original_transfer': 'reuse_selected_ema'}
            plan = {'schema': lr.SCHEMA, 'seed': 17, 'fold': 0, 'common_step_policy': policy,
                    'original_recipe': {'config': config}, 'recipe': variant, 'recipe_sha256': json_hash(variant),
                    'primary_screen': None, 'baseline': baseline_binding, 'recipient': bind(recipient),
                    'mapped_state_sha256': 'fixture'}
            plan_path = slot / 'plan.json'
            write_json(plan_path, plan)
            plan_sha = base.sha256(plan_path)
            directory = lr.checkpoint_root(runtime, 17)
            directory.mkdir(parents=True)
            checkpoint = directory / lr.NATIVE_CHECKPOINT
            checkpoint.touch()
            metadata = {'identity': {**full['identity'], 'arm': lr.ARM, 'regime': base.REGIME,
                                     'seed': 17, 'stage': 'formal'},
                        'recipe_sha256': json_hash(variant), 'best_step': 4500, 'best_validation_loss': .5,
                        'model_state_sha256': 'fixture', 'checkpoint_sha256': base.sha256(checkpoint)}
            history = [{'step': 4500, 'validation_loss': .5, 'early_stop': True}]
            write_json(directory / 'native-selection.json', {**metadata, 'monitor': MONITOR, 'stop_step': 4500,
                       'stop_reason': 'validation_early_stop', 'history_sha256': json_hash(history)})
            write_json(directory / 'run.json', {**metadata, 'recipe': variant, 'monitor': MONITOR,
                       'history': history, 'status': 'incomplete', 'optimizer_steps': 6000})
            write_json(directory / 'fit.json', {'plan_sha256': plan_sha, 'available_modes': ['native'],
                       'initialization': {'mapped_state_sha256': 'fixture', 'mapped_dose_rows_trainable': True},
                       'files': {path.name: base.sha256(path) for path in directory.iterdir()}})
            self.write_prediction(slot / 'predictions/native',
                                  {'checkpoint': bind(checkpoint), 'metadata': metadata}, request,
                                  {'plan_sha256': plan_sha, 'arm': lr.ARM, 'mode': 'native'}, [[1, -1]])

            self.cli('score', root, '--runtime', runtime, '--plan-sha256', plan_sha)
            report = base.read_json(slot / 'scores/scores.json')
            self.assertEqual(report['comparisons']['common']['status'], 'incomplete')
            self.assertEqual(report['comparisons']['native']['paired'][0]['mean_ap100']['gain'], .5)
            self.assertEqual(report['comparisons']['native']['paired'][0]['mean_ap100']['recovery'], 1)
            self.cli('decision', root)
            self.assertEqual(base.read_json(root / 'decision.json')['status'], 'incomplete')

    def write_prediction(self, destination, item, request, lineage, values):
        cases = base.case_order(base.case_rows(request))
        write_json(destination / 'cases.json', cases)
        prediction = destination / 'prediction.npy'
        base.np.save(prediction, base.np.array(values, dtype='float32'))
        write_json(destination / 'prediction.json', {**item, 'schema': base.PREDICTION_SCHEMA,
                   'lineage': lineage, 'request_sha256': json_hash(request), 'case_count': 1,
                   'group_count': 1, 'ordered_case_group_sha256': json_hash(cases),
                   'cases_sha256': base.sha256(destination / 'cases.json'),
                   'prediction_sha256': base.sha256(prediction)})

    def cli(self, phase, root, *args):
        result = subprocess.run([sys.executable, lr.__file__, phase, str(root), *map(str, args)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
