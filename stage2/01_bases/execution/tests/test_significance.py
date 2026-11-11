"""Published truth and partition access contracts, independent numeric examples."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.significance import ranked_indices


class TruthContract(unittest.TestCase):
    def test_directional_lists_use_padj_effect_then_integer_index(self):
        up, down = ranked_indices(
            np.array([.01, .01, .01, .03, .05, np.nan]),
            np.array([2., 3., 3., -4., -8., 9.]),
        )
        self.assertEqual(up, [1, 2, 0])
        self.assertEqual(down, [3])


if __name__ == '__main__':
    unittest.main()

class ContrastMembershipContract(unittest.TestCase):
    def test_source_forms_join_parent_dose_with_exact_plate_controls(self):
        import pandas as pd
        from bases.significance import contrast_plan
        native = pd.DataFrame({
            'plate': ['P1'] * 4, 'is_control': [False, False, True, True],
            'is_blank': [False] * 4, 'compound_name': ['salt', 'base', 'DMSO', 'Media'],
            'vehicle': ['DMSO', 'DMSO', 'DMSO', 'Media'], 'time_h': [24] * 4,
            'drug_key': ['parent', 'parent', '', ''], 'dose_um': [1, 1, 0, 0],
        }, index=['well1', 'well2', 'dmso', 'media'])
        treated = pd.DataFrame({
            'sample_id': ['s1', 's2'], 'source_profile_id': ['well1', 'well2'],
            'source_row': [0, 1], 'drug_parent': ['parent', 'parent'], 'dose_um': [1., 1.],
            'matched_control_rows': ['[2]', '[2]'],
            'matched_control_ids': ['["dmso"]', '["dmso"]'],
        })
        groups, samples = contrast_plan(treated, native)
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups.iloc[0].n_treated, 2)
        self.assertEqual(groups.iloc[0].n_controls, 1)
        self.assertEqual(samples.source_profile_id.tolist(), ['well1', 'well2', 'dmso'])
        native.loc['well2', 'time_h'] = 48
        with self.assertRaisesRegex(ValueError, 'incompatible'):
            contrast_plan(treated, native)
        native.loc['well2', 'time_h'] = 24
        treated.loc[1, 'matched_control_rows'] = '[3]'
        with self.assertRaisesRegex(ValueError, 'DMSO'):
            contrast_plan(treated, native)

class PartitionTruthContract(unittest.TestCase):
    def test_release_distinguishes_empty_untestable_and_rejects_split_group(self):
        import tempfile
        import pandas as pd
        from bases.significance import group_key, publish_statistics, PartitionTruth
        groups = pd.DataFrame([
            {'group_id': group_key(parent, 1), 'drug_parent': parent, 'dose_um': 1,
             'status': status, 'reason': reason, 'sample_ids_sha256': digest}
            for parent, status, reason, digest in [
                ('A', 'eligible', '', __import__('hashlib').sha256(b'["a1","a2"]').hexdigest()),
                ('B', 'eligible', '', __import__('hashlib').sha256(b'["b1"]').hexdigest()),
                ('C', 'untestable', 'insufficient_replication', __import__('hashlib').sha256(b'["c1"]').hexdigest()),
            ]])
        samples = pd.DataFrame([
            {'group_id': group_key(parent, 1), 'sample_id': sample,
             'source_profile_id': 'native-' + sample, 'condition': 'treated'}
            for parent, sample in [('A', 'a1'), ('A', 'a2'), ('B', 'b1'), ('C', 'c1')]])
        stats = {'pvalue': np.array([[.001,.5,.5], [.5,.5,.5], [np.nan]*3]),
                 'padj': np.array([[.003,.5,.5], [.5,.5,.5], [np.nan]*3]),
                 'log2_fold_change': np.array([[2,1,-1], [1,1,-1], [np.nan]*3])}
        genes = pd.DataFrame({'gene_index':[0,1,2], 'hgnc_id':['HGNC:1','HGNC:2','HGNC:3']})
        with tempfile.TemporaryDirectory() as directory:
            release = Path(directory) / 'truth'
            publish_statistics(release, genes, groups, samples, stats, {}, fixture=True)
            assignment = pd.DataFrame({'sample_id':['a1','a2','b1','c1'],
                                       'role':['test','test','validation','test']})
            labels = PartitionTruth(release, assignment, 'test', allow_fixture=True)
            self.assertEqual(labels.truth('A', 1, ['a2','a1']).tolist(), [1,0,0])
            self.assertIsNone(labels.truth('C', 1, ['c1']))
            self.assertEqual(labels.ranked('A', 1), ([0], []))
            with self.assertRaisesRegex(ValueError, 'partition'):
                labels.truth('B', 1, ['b1'])
            empty = PartitionTruth(release, assignment, 'validation', allow_fixture=True)
            self.assertEqual(empty.truth('B', 1, ['b1']).tolist(), [0,0,0])
            assignment.loc[1, 'role'] = 'train'
            with self.assertRaisesRegex(ValueError, 'crosses'):
                PartitionTruth(release, assignment, 'test', allow_fixture=True)
