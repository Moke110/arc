"""Final Step 06 attribution and paired publication contracts."""
from pathlib import Path
import tempfile
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
import publication
from publication import diagnosis


class PublicationTest(unittest.TestCase):
    def test_attribution_requires_all_gates_on_the_primary_comparison(self):
        screen = {'gains': {'17': '0.010', '42': '0.001', '2026': '0.002'}}
        metric = {'gain': 0.02, 'cm_gap': 0.04, 'recovery': 0.5,
                  'paired_seeds': [17, 42, 2026], 'matched_cm_seeds': [17, 42, 2026]}
        arm = {'primary_mode': 'common',
               'comparisons': {'common': {'pooled': {'mean_ap100': metric}}}}
        self.assertTrue(diagnosis(screen, arm)['major_contribution'])
        for change, reason in [
            ({'gain': 0}, 'nonpositive_confirmation_effect'),
            ({'cm_gap': 0, 'recovery': None}, 'nonpositive_cm_gap'),
            ({'recovery': None, 'cm_gap': None, 'matched_cm_seeds': []}, 'matching_cm_unavailable'),
            ({'recovery': 0.4999}, 'below_half_gap'),
            ({'paired_seeds': [17, 42]}, 'incomplete_confirmation'),
        ]:
            with self.subTest(reason=reason):
                metric_before = metric.copy()
                metric.update(change)
                result = diagnosis(screen, arm)
                self.assertFalse(result['major_contribution'])
                self.assertEqual(result['attribution_status'], reason)
                metric.clear()
                metric.update(metric_before)

    def test_unpromoted_arms_retain_their_distinct_terminal_reasons(self):
        for gains, status in [
            ({'17': '0.009999'}, 'below_threshold'),
            ({'17': '0.010', '42': '0.001', '2026': '0'}, 'inconsistent'),
            ({'17': '0.010', '42': '0.001', '2026': None}, 'incomplete'),
        ]:
            with self.subTest(status=status):
                result = diagnosis({'gains': gains}, None)
                self.assertEqual(result['attribution_status'], status)
                self.assertFalse(result['major_contribution'])
                self.assertIsNone(result['confirmation'])

    def test_record_and_figures_replay_and_detect_changed_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output, record = root / 'output', root / 'record'
            panel = output / 'source-audit' / 'exposure.csv'
            exposure = [{'source_release': source, 'tahoe_role': role, 'tenv_role': target,
                         'shared_parent_count': int(source == 'historical' and role == 'train'
                                                    and target == 'validation')}
                        for source in ('historical', 'globally_cold')
                        for role in ('train', 'validation')
                        for target in ('train', 'validation', 'test')]
            publication.write_csv(panel, exposure)
            original = [{'dose_um': 'pooled', 'group_count': 1, 'positive_ap100': .2,
                         'negative_ap100': .4, 'mean_ap100': .3}]
            treatment = [{**original[0], 'positive_ap100': .22, 'mean_ap100': .31}]
            cm = [{**original[0], 'positive_ap100': .24, 'mean_ap100': .32}]
            pairs = list(publication.pair_rows(0, 17, 'dose_reset', 'loss-selected',
                                               {'status': 'completed'}, original, treatment, cm,
                                               (6500, 6000, 6500)))
            summary = {'replication': {'status': 'skipped'}}
            with patch.object(publication, 'collect', return_value=(summary, pairs, {})):
                publication.write_record(root, output, record, False)
                self.assertEqual(publication.write_record(root, output, record, True), summary)
                rows = publication.rows(record / 'per_run.csv')
                self.assertEqual({row['metric'] for row in rows},
                                 {'mean_ap100', 'positive_ap100', 'negative_ap100'})
                pooled = next(row for row in rows if row['metric'] == 'mean_ap100')
                self.assertEqual((pooled['group_count'], pooled['gain'], pooled['recovery']),
                                 ('1', '0.01', '0.5'))
                figure = output / publication.FIGURES[0]
                self.assertTrue(figure.read_bytes().startswith(b'\x89PNG'))
                figure.write_bytes(b'changed PNG')
                with self.assertRaisesRegex(ValueError, publication.FIGURES[0]):
                    publication.write_record(root, output, record, True)


if __name__ == '__main__':
    unittest.main()
