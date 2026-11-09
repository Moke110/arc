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


def promoted_fixture(output):
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
    upstream = {'status': 'ready', 'seeds': list(SEEDS), 'promoted_arms': ['learning_rate'], 'recipient': target,
                'comparison': {'prediction': {}, 'score': 'directional mean AP@100',
                               'common_steps': {str(seed): 3 for seed in SEEDS}},
                'screens': {'dose_reset': {'status': 'below_threshold'}, 'learning_rate': {'status': 'promoted'}}}
    for seed in SEEDS:
        write_json(replication.confirmation.source_path(output, seed) / 'source.json', {'seed': seed})
    original = {'learning_rate': 0.0002}
    cm = {'config': original, 'metadata': {'best_step': 2, 'identity': {'fold': 1}}}
    with (patch.object(replication.confirmation, 'load_plan', return_value=upstream),
          patch.object(replication.confirmation, 'load_source', return_value={'checkpoint_sha256': 'tahoe'}),
          patch.object(replication, 'historical_scoring', return_value={'step04_checkpoints': bind(manifest)}),
          patch.object(replication, 'matched_cm', return_value=cm),
          patch.object(replication.baseline, 'verify_inference', return_value=cases)):
        digest = replication.freeze(output, output, 'confirmation-sha')['sha256']
        frozen = replication.load_plan(output, output, digest)
        for seed in SEEDS:
            item = frozen['seeds'][str(seed)]
            assert item['arms']['learning_rate']['config'] == {'learning_rate': .0001}
            assert item['arms']['Transfer']['config'] == original
            slot = replication.seed_slot(output, seed)
            chosen, fits = {}, {}
            for arm in ('CM', 'Transfer', 'learning_rate'):
                if arm != 'CM':
                    fit = replication.fit_binding(output, seed, arm) / 'fit.json'
                    write_json(fit, {'fixture': arm})
                    fits[arm] = {'status': 'completed', 'reason': 'max_steps', 'checkpoints': {}}
                else:
                    fit = manifest
                for mode in (('native',) if arm == 'CM' else ('native', 'common')):
                    key = f'{arm}/{mode}'
                    selection = {'checkpoint': bind(fit), 'metadata': {}, 'fit': bind(fit)}
                    chosen[key] = selection
                    dest = slot / 'predictions' / key
                    dest.mkdir(parents=True)
                    prediction = -truth if (arm == 'Transfer') != (seed == 2026) else truth
                    np.save(dest / 'prediction.npy', prediction)
                    order = replication.baseline.case_order(cases)
                    write_json(dest / 'cases.json', order)
                    write_json(dest / 'prediction.json', {
                        'schema': replication.baseline.PREDICTION_SCHEMA,
                        'lineage': replication.prediction_lineage(digest, seed, key, selection),
                        'checkpoint': selection['checkpoint'], 'metadata': {},
                        'request_sha256': json_hash(inference), 'case_count': 1, 'group_count': 1,
                        'ordered_case_group_sha256': json_hash(order),
                        'cases_sha256': bind(dest / 'cases.json')['sha256'],
                        'prediction_sha256': bind(dest / 'prediction.npy')['sha256'],
                    })
            with (patch.object(replication, 'comparisons_for', return_value=(chosen, fits)),
                  patch.object(replication.baseline, 'case_rows', return_value=cases)):
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
