"""Bounded fold 1 replication contracts."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
import replication


class ReplicationTest(unittest.TestCase):
    def test_skipped_confirmation_needs_no_fold1_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            plan = {'status': 'skipped', 'seeds': [], 'promoted_arms': [],
                    'screens': {'dose_reset': {'status': 'below_threshold'},
                                'learning_rate': {'status': 'below_threshold'}}}
            with patch.object(replication.confirmation, 'load_plan', return_value=plan):
                binding = replication.freeze(output, output, 'confirmation-sha')
                frozen = replication.load_plan(output, output, binding['sha256'])
                report = replication.finalize(output, frozen, binding['sha256'], {})
            self.assertEqual(report['status'], 'skipped')
            self.assertEqual(report['screens'], plan['screens'])
            self.assertFalse((output / 'checkpoints').exists())
            self.assertFalse((output / 'replication/fold-1-seed-17').exists())

    def test_positive_equal_step_pair_without_matched_cm_has_no_recovery(self):
        left = {'dose_um': 'pooled', 'group_count': 2, 'positive_ap100': .2,
                'negative_ap100': .4, 'mean_ap100': .3}
        right = {'dose_um': 'pooled', 'group_count': 2, 'positive_ap100': .4,
                 'negative_ap100': .4, 'mean_ap100': .4}
        summaries = {'Transfer/common': [left], 'learning_rate/common': [right],
                     'CM/native': [{**right, 'mean_ap100': .5}]}
        row = replication.compare(summaries, 'learning_rate', 'common', False)
        metric = row['contrasts'][0]['mean_ap100']
        self.assertEqual(metric['gain'], .1)
        self.assertIsNone(metric['cm_gap'])
        self.assertIsNone(metric['recovery'])
        self.assertIn('matching CM', row['cm_reason'])

    def test_promoted_plan_scores_all_pairs_and_keeps_cm_modes_separate(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            report = promoted_fixture(output)
            arm = report['arms']['learning_rate']
            self.assertEqual(report['status'], 'completed')
            self.assertTrue(arm['replicated'])
            pooled = arm['comparisons']['common']['pooled']['mean_ap100']
            self.assertEqual(pooled['paired_seeds'], [17, 42, 2026])
            self.assertEqual(pooled['matched_cm_seeds'], [])
            self.assertAlmostEqual(pooled['gain'], 2 / 9)
            negative = arm['comparisons']['common']['per_seed']['2026']['contrasts'][0]['mean_ap100']
            self.assertAlmostEqual(negative['gain'], -2 / 3)
            self.assertIsNone(pooled['recovery'])
            native = arm['comparisons']['native']['pooled']['mean_ap100']
            self.assertEqual(native['matched_cm_seeds'], [17, 42, 2026])
            self.assertEqual(native['recovery'], 1)

    def test_failed_source_keeps_surviving_seed_pairs_but_cannot_replicate(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = promoted_fixture(Path(tmp), source_unavailable={
                '2026': ('failed', 'Tahoe source fit exceeded the resource allowance')})
            arm = report['arms']['learning_rate']
            self.assertEqual(report['status'], 'incomplete')
            self.assertIsNone(arm['replicated'])
            pooled = arm['comparisons']['common']['pooled']['mean_ap100']
            self.assertEqual(pooled['paired_seeds'], [17, 42])
            self.assertEqual(pooled['planned_seed_count'], 3)
            self.assertEqual(report['seed_dispositions']['2026'], {
                'status': 'failed', 'reason': 'Tahoe source fit exceeded the resource allowance'})

    def test_failed_target_arm_keeps_other_intervention_pair(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = promoted_fixture(Path(tmp), arm_unavailable={
                'dose_reset': ('failed', 'target fit interrupted before sealing')})
            lr = report['arms']['learning_rate']
            reset = report['arms']['dose_reset']
            self.assertEqual(report['status'], 'incomplete')
            self.assertTrue(lr['replicated'])
            self.assertEqual(lr['comparisons']['common']['pooled']['mean_ap100']['paired_seeds'],
                             [17, 42, 2026])
            self.assertIsNone(reset['replicated'])
            self.assertEqual(reset['comparisons']['native']['pooled']['mean_ap100']['paired_seeds'],
                             [42, 2026])
            pair = reset['comparisons']['native']['per_seed']['17']
            self.assertIn('target fit interrupted before sealing', pair['reason'])


def promoted_fixture(output, source_unavailable=None, arm_unavailable=None):
    from artifacts import SEEDS, write_csv, write_json
    from formal_matrix import bind
    from training_inputs import json_hash

    np = replication.np
    cases = [{'case_id': 'well-1', 'group_id': 'drug1-dose1', 'drug_key': 'drug1', 'dose_um': '1.0'}]
    truth = np.array([[1., 0., -1.]], dtype='float32')
    groups = [{key: cases[0][key] for key in ('group_id', 'drug_key', 'dose_um')} | {'n_wells': 1}]
    prefix = output / 'fold_1/evaluator'
    write_csv(prefix / 'groups.csv', groups)
    np.save(prefix / 'truth.npy', truth)
    for name, label in zip(('positive', 'negative'), replication.baseline.truth_labels(truth[0])):
        np.save(prefix / f'{name}.npy', label[None, :])
    manifest = output / 'manifest.json'
    write_json(manifest, {'files': {str(path.relative_to(output)): bind(path)['sha256']
                                   for path in sorted(prefix.iterdir())}})
    request = {'input_root': str(output), 'identity': {'source_input_sha256': bind(manifest)['sha256']}}
    inference = {'n_genes': 3}
    target = {'request': request, 'inference': inference, 'producer_request': bind(manifest),
              'fixed_config': bind(manifest), 'roles': {'validation': {'group_order_sha256': 'val'}}}
    arms = ['learning_rate', 'dose_reset'] if arm_unavailable else ['learning_rate']
    upstream = {'status': 'ready', 'seeds': list(SEEDS), 'promoted_arms': arms, 'recipient': target,
                'comparison': {'prediction': {}, 'score': 'directional mean AP@100',
                               'common_steps': {str(seed): 3 for seed in SEEDS}},
                'screens': {'dose_reset': {'status': 'below_threshold'}, 'learning_rate': {'status': 'promoted'}}}
    for seed in SEEDS:
        if str(seed) not in (source_unavailable or {}):
            write_json(replication.confirmation.source_path(output, seed) / 'source.json', {'seed': seed})
    original = {'learning_rate': 0.0002}
    cm = {'config': original, 'metadata': {'best_step': 2, 'identity': {'fold': 1}},
          **{key: bind(manifest) for key in ('history', 'replay', 'recipe', 'request', 'checkpoint')}}
    with (patch.object(replication.confirmation, 'load_plan', return_value=upstream),
          patch.object(replication.confirmation, 'load_source', return_value={'checkpoint_sha256': 'tahoe'}),
          patch.object(replication, 'historical_scoring', return_value={'step04_checkpoints': bind(manifest)}),
          patch.object(replication, 'matched_cm', return_value=cm),
          patch.object(replication.baseline, 'verify_inference', return_value=cases)):
        digest = (replication.freeze(output, output, 'confirmation-sha', source_unavailable) if source_unavailable
                  else replication.freeze(output, output, 'confirmation-sha'))['sha256']
        frozen = replication.load_plan(output, output, digest)
        for seed in SEEDS:
            if str(seed) in (source_unavailable or {}):
                continue
            item = frozen['seeds'][str(seed)]
            assert item['arms']['learning_rate']['config'] == {'learning_rate': .0001}
            assert item['arms']['Transfer']['config'] == original
            slot = replication.seed_slot(output, seed)
            fits = {}
            available = ['CM', 'Transfer', 'learning_rate']
            if arm_unavailable and seed != 17:
                available.append('dose_reset')
            for arm in available:
                if arm != 'CM':
                    fit = replication.fit_binding(output, seed, arm) / 'fit.json'
                    write_json(fit, {'fixture': arm})
                    fits[arm] = {'status': 'completed', 'reason': 'max_steps', 'checkpoints': {},
                                 'bridge': bind(manifest)}
                else:
                    fit = manifest
                for mode in (('native',) if arm in ('CM', 'dose_reset') else ('native', 'common')):
                    key = f'{arm}/{mode}'
                    metadata = cm['metadata'] if arm == 'CM' else {'best_step': 3}
                    selection = {'checkpoint': bind(fit), 'metadata': metadata, 'fit': bind(fit)}
                    if arm != 'CM':
                        fits[arm]['checkpoints'][mode] = selection
                    dest = slot / 'predictions' / key
                    dest.mkdir(parents=True)
                    prediction = -truth if (arm == 'Transfer') != (seed == 2026) else truth
                    np.save(dest / 'prediction.npy', prediction)
                    order = replication.baseline.case_order(cases)
                    write_json(dest / 'cases.json', order)
                    write_json(dest / 'prediction.json', {
                        'schema': replication.baseline.PREDICTION_SCHEMA,
                        'lineage': replication.prediction_lineage(digest, seed, key, selection),
                        'checkpoint': selection['checkpoint'], 'metadata': metadata,
                        'request_sha256': json_hash(inference), 'case_count': 1, 'group_count': 1,
                        'ordered_case_group_sha256': json_hash(order),
                        'cases_sha256': bind(dest / 'cases.json')['sha256'],
                        'prediction_sha256': bind(dest / 'prediction.npy')['sha256'],
                    })
            with (patch.object(replication, 'fit_selected', side_effect=lambda *args: fits[args[-1]]),
                  patch.object(replication.baseline, 'predict_checkpoint'),
                  patch.object(replication, 'use_runtime'),
                  patch.object(replication.baseline, 'case_rows', return_value=cases)):
                if arm_unavailable and seed == 17:
                    replication.predict(output, output, frozen, digest, seed, arm_unavailable)
                prediction = slot / 'predictions/Transfer/common/prediction.npy'
                saved = prediction.read_bytes()
                prediction.write_bytes(b'wrong prediction')
                # A bad prediction fails before opening the intentionally missing truth.
                observed = prefix / 'truth.npy'
                moved = prefix / 'temporarily-missing.npy'
                observed.rename(moved)
                try:
                    replication.score(output, frozen, digest, seed)
                except ValueError as error:
                    assert 'prediction.npy' in str(error), str(error)
                else:
                    raise AssertionError('tampered prediction was accepted')
                moved.rename(observed)
                prediction.write_bytes(saved)
                replication.score(output, frozen, digest, seed)
        return replication.finalize(output, frozen, digest, {})


if __name__ == '__main__':
    unittest.main()
