import shutil
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from native_log_dose_fixture import source_release
import paired_reference as reference
import paired_target as target
from artifacts import SEEDS, read_json, write_csv, write_json
from formal_matrix import bind
from paired_evidence import verify_matrix
from paired_log_source import source_inputs, source_slot
from training_inputs import json_hash
from training_policy import MONITOR


def fixture_fit(directory, digest, identity, config, groups, losses, **extra):
    directory.mkdir(parents=True)
    checkpoint = directory / 'best.pt'
    checkpoint.write_bytes(b'fixture EMA')
    step, loss = min(enumerate(losses, 1), key=lambda item: item[1])
    metadata = {'identity': identity, 'recipe_sha256': json_hash(config),
                'checkpoint_sha256': bind(checkpoint)['sha256'], 'best_step': step,
                'best_validation_loss': loss, 'model_state_sha256': 'fixture-state'}
    write_json(directory / 'run.json', {
        **metadata, **extra, 'status': 'completed', 'monitor': MONITOR, 'recipe': config,
        'optimizer_steps': len(losses), 'stop_reason': 'max_steps',
        'history': [{'step': i, 'validation_loss': value, 'early_stop': False}
                    for i, value in enumerate(losses, 1)],
        'reload': {'validation_loss': loss, 'prediction_max_abs_error': 0},
        'validation_batching': {'group_ids': groups, 'group_ids_sha256': json_hash(groups)},
    })
    history = bind(directory / 'run.json')
    write_json(directory / 'fit.json', {'plan_sha256': digest, 'run_sha256': history['sha256']})
    return {'metadata': metadata, 'checkpoint': bind(checkpoint), 'history': history,
            'fit_seal': bind(directory / 'fit.json')}


def fixture_score(root, plan, digest, arm, selected, predictions):
    destination = root / ('Transfer' if arm == 'Transfer' else 'prediction')
    destination.mkdir(parents=True)
    np.save(destination / 'prediction.npy', np.asarray(predictions, dtype='float32'))
    cases = reference.base.case_order(reference.base.case_rows(plan['request']))
    write_json(destination / 'cases.json', cases)
    lineage = {'plan_sha256': digest, 'protocol_sha256': plan['protocol_sha256'], 'arm': arm}
    write_json(destination / 'prediction.json', {
        'schema': reference.base.PREDICTION_SCHEMA, 'lineage': lineage, 'checkpoint': selected['checkpoint'],
        'metadata': selected['metadata'], 'request_sha256': json_hash(plan['request']),
        'case_count': len(cases), 'group_count': plan['group_count'],
        'ordered_case_group_sha256': json_hash(cases), 'cases_sha256': bind(destination / 'cases.json')['sha256'],
        'prediction_sha256': bind(destination / 'prediction.npy')['sha256'],
    })
    return reference.score_paired(root, plan, digest, arm=arm, selected=selected,
                                  prediction_dir=destination, lineage=lineage)


def fixture_matrix(root):
    shared = root / 'inputs'
    shared.mkdir()
    files = {key: shared / key for key in reference.base.INFERENCE_FILES}
    for file in files.values():
        file.touch()
    doses = ['0.2', '1.0', '10.0']
    groups = ['a', 'b', 'c']
    cases = [{'case_id': f'case-{group}', 'group_id': group, 'dose_um': dose, 'drug_key': group}
             for group, dose in zip(groups, doses)]
    for role in ('train', 'validation', 'test'):
        write_csv(files[f'{role}_cases'], cases)
    truth = shared / 'truth.npy'
    np.save(truth, np.array([[1, -1, 0], [1, 1, 1], [1, -1, 0]], dtype='float32'))
    identity = {'dataset': 'tenv', 'fold': 0, 'gene_sha256': bind(files['genes'])['sha256']}
    request = {'identity': identity, 'n_genes': 3, 'target_kind': 'delta',
               'predictor_files': {k: str(v) for k, v in files.items()},
               'file_sha256': {k: bind(v)['sha256'] for k, v in files.items()}}
    frozen = {'fold': 0, 'regime': reference.REGIME, 'request': reference.base.inference_request(request),
              'case_count': 3, 'group_count': 3, 'dose_group_counts': dict.fromkeys(doses, 1),
              'ordered_case_group_sha256': json_hash(reference.base.case_order(cases)), 'truth': bind(truth)}
    write_json(shared / 'request.json', request)
    write_json(shared / 'recipe.json', {})
    matrix = root / 'output' / 'paired-transfer'
    for seed in SEEDS:
        upstream = matrix / 'reference' / f'seed-{seed}'
        out = matrix / 'log' / f'seed-{seed}'
        config = {'seed': seed, 'stage': 'formal', 'max_steps': 3}
        original_source = source_release(SimpleNamespace(np=np), root / f'tahoe-{seed}', config)
        original_identity = {**identity, 'arm': 'Transfer', 'regime': reference.REGIME,
                             'seed': seed, 'stage': 'formal'}
        original = fixture_fit(upstream / 'fit', 'step05-plan', original_identity, config, groups, [.7, .5, .6])
        plan = {**frozen, 'schema': reference.SCHEMA, 'seed': seed, 'config': config,
                'source': original_source, 'protocol_sha256': json_hash(frozen),
                'recipient_request': bind(shared / 'request.json'), 'target_recipe': bind(shared / 'recipe.json'),
                'checkpoints': {'Transfer': original}}
        write_json(upstream / 'plan.json', plan)
        reference_digest = bind(upstream / 'plan.json')['sha256']
        original_score = fixture_score(upstream, plan, reference_digest, 'Transfer', original,
                                       [[-1, 1, 0], [0, 0, 0], [-1, 1, 0]])
        _, _, _, source_protocol = source_inputs(plan, reference_digest)
        fixture_fit(source_slot(root, seed) / 'fit', json_hash(source_protocol),
                    source_protocol['identity'], config, ['validation'], [.8, .6, .6])
        target_plan_file = target.freeze(out, 'log', seed, upstream, reference_digest, runtime=root)
        digest = target_plan_file['sha256']
        target_plan = read_json(target_plan_file['path'])
        selection = fixture_fit(Path(target_plan['checkpoint_slot']) / 'fit', digest,
                                target.fit_identity(target_plan, plan), config, groups, [.7, .8, .9],
                                source=target_plan['source'], training_initial_state={
                                    'optimizer_state_entries': 0, 'scheduler_start_step': 0, 'ema_start_step': 0})
        scores = fixture_score(out, plan, digest, 'log', selection, [[1, -1, 0], [0, 0, 0], [0, 0, 0]])
        write_json(out / 'paired.json', {
            'schema': target.PAIR_SCHEMA, 'arm': 'log', 'seed': seed, 'plan_sha256': digest,
            'protocol_sha256': plan['protocol_sha256'], 'scores': scores, 'reference': original_score,
            'per_dose': target.paired_differences(read_json(scores['path']), read_json(original_score['path'])),
        })
    return matrix


class PairedEvidence(unittest.TestCase):
    def test_matrix_keeps_loss_selections_and_undefined_support_and_requires_all_seeds(self):
        with tempfile.TemporaryDirectory() as temp:
            root = fixture_matrix(Path(temp))
            result = verify_matrix(root, 'log', runtime=Path(temp))
            self.assertEqual([cell['seed'] for cell in result['cells']], [17, 42, 2026])
            for cell in result['cells']:
                self.assertEqual(cell['source']['metadata']['best_step'], 2)
                self.assertEqual(cell['target']['metadata']['best_step'], 1)
                self.assertEqual(cell['optimizer_steps'], 3)
                self.assertAlmostEqual(cell['per_dose'][0]['pcc_difference'], 2)
                self.assertEqual(cell['per_dose'][1]['pcc_undefined_reasons'], {'constant_truth': 1})
                self.assertEqual(cell['per_dose'][2]['pcc_undefined_reasons'], {'constant_prediction': 1})
                self.assertIsNone(cell['per_dose'][2]['pcc_difference'])
            (root / 'log' / 'seed-2026' / 'paired.json').unlink()
            with self.assertRaisesRegex(FileNotFoundError, 'seed-2026'):
                verify_matrix(root, 'log', runtime=Path(temp))

    def test_matrix_rejects_a_seal_substituted_from_another_seed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = fixture_matrix(Path(temp))
            shutil.copyfile(root / 'log' / 'seed-17' / 'paired.json', root / 'log' / 'seed-42' / 'paired.json')
            with self.assertRaisesRegex(ValueError, 'paired score lineage'):
                verify_matrix(root, 'log', runtime=Path(temp))

    def test_matrix_rejects_scores_from_before_a_prediction_was_rebound(self):
        with tempfile.TemporaryDirectory() as temp:
            root = fixture_matrix(Path(temp))
            cell = root / 'log/seed-17'
            values = cell / 'prediction/prediction.npy'
            np.save(values, np.array([[-1, 1, 0], [0, 0, 0], [0, 0, 0]], dtype='float32'))
            prediction = values.with_name('prediction.json')
            write_json(prediction, {**read_json(prediction), 'prediction_sha256': bind(values)['sha256']})
            report = cell / 'scores/scores.json'
            write_json(report, {**read_json(report), 'prediction': bind(prediction)})
            pair = cell / 'paired.json'
            write_json(pair, {**read_json(pair), 'scores': bind(report)})
            with self.assertRaisesRegex(ValueError, 'group scores differ'):
                verify_matrix(root, 'log', runtime=Path(temp))
