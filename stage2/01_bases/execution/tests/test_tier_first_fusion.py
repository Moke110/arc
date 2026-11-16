"""ADR0073 behavior through the existing public disease fusion boundary."""
import sys
import unittest
import tempfile
import hashlib
import json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api
import build_disease_fusion as driver
from test_disease_fusion import write_fixture, context_for
from bases import disease_deg

class TierFirstFusionTests(unittest.TestCase):
    def test_actual_higher_tier_deg_support_precedes_better_numeric_sum(self):
        result = disease_deg.fuse_disease([
            {'dataset_id': 'strong', 'tier': 1, 'n_deg': 2,
             'rank': np.array([2, 1, 3, 3, 3])},
            {'dataset_id': 'weak', 'tier': 2, 'n_deg': 3,
             'rank': np.array([4, 2, 1, 3, 4])},
        ], axis_size=5)
        self.assertGreater(result['sum'][0], result['sum'][2])
        self.assertLess(result['rank'][0], result['rank'][2])
        self.assertEqual(result['source_tier_highest'].tolist(), [1, 1, 2, 2, None])

    def test_cutoff_retains_complete_union_and_supporting_source_evidence(self):
        result = disease_deg.fuse_disease([
            {'dataset_id': 'one', 'tier': 3, 'n_deg': 1001,
             'rank': np.arange(1, 1003)},
        ], axis_size=1002)
        self.assertEqual(result['deg1000'], 1000)
        full = disease_deg.union_rows(result)
        kept = disease_deg.deg1000_rows(result)
        self.assertEqual(len(full), 1001)
        self.assertEqual(len(kept), 1000)
        self.assertEqual(full.iloc[-1].disease_rank, 1001)
        self.assertEqual(full.iloc[-1].source_tier_highest, 3)
        self.assertEqual(full.iloc[-1].supporting_sources, '["one"]')
        self.assertEqual(full.iloc[-1].disease_rank_sum, 1.0)

    def test_published_reader_exposes_full_union_with_tier_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root, sources = write_fixture(directory)
            context = context_for(root, sources)
            driver.publish(context)
            reader = disease_deg.DiseaseFusion(context.output)
            full = reader.union('nhf')
            self.assertEqual(full.gene_index.tolist(), [2, 1])
            self.assertEqual(full.source_tier_highest.tolist(), [1, 1])
            self.assertEqual(full.supporting_sources.tolist(), ['["nhf_a"]', '["nhf_a"]'])
            self.assertEqual(driver.replay(context)['result'], 'PASS')

    def test_missing_or_invalid_source_tier_fails_closed(self):
        for tier in [None, 0, 4, True, 1.0, '1']:
            with self.subTest(tier=tier), self.assertRaises(ValueError):
                disease_deg.fuse_disease([
                    {'dataset_id': 'bad', 'tier': tier, 'n_deg': 1,
                     'rank': np.array([1, 2])}], axis_size=2)

    def test_gate_and_cutoff_boundaries_retain_untruncated_ranks(self):
        for count in [0, 99, 100, 999, 1000, 1001]:
            with self.subTest(count=count):
                rank = np.full(1002, count + 1)
                rank[:count] = np.arange(1, count + 1)
                result = disease_deg.fuse_disease([
                    {'dataset_id': 'source', 'tier': 1, 'n_deg': count,
                     'rank': rank}], axis_size=1002)
                self.assertEqual(result['eligible'], count >= 100)
                self.assertEqual(result['deg1000'], min(count, 1000))
                full = disease_deg.union_rows(result)
                self.assertEqual(full.disease_rank.tolist(), list(range(1, count + 1)))
                self.assertEqual(len(disease_deg.deg1000_rows(result)), min(count, 1000))

    def test_noninformative_source_cannot_confer_tier_and_diseases_are_independent(self):
        weak = {'dataset_id': 'weak', 'tier': 3, 'n_deg': 1, 'rank': np.array([1, 2])}
        excluded = {'dataset_id': 'excluded', 'tier': 1, 'n_deg': 1,
                    'rank': np.array([1, 2]), 'informative': False}
        first = disease_deg.fuse_disease([weak, excluded], axis_size=2)
        second = disease_deg.fuse_disease([dict(weak, tier=2)], axis_size=2)
        self.assertEqual(first['source_tier_highest'].tolist(), [3, None])
        self.assertEqual(second['source_tier_highest'].tolist(), [2, None])
        self.assertEqual(first['supporting_sources'][0], ['weak'])

    def test_reader_rejects_forged_tier_even_when_file_digest_is_updated(self):
        with tempfile.TemporaryDirectory() as directory:
            root, sources = write_fixture(directory)
            context = context_for(root, sources)
            driver.publish(context)
            reader = disease_deg.DiseaseFusion(context.output)
            path = context.output / reader.record('nhf')['union_file']
            frame = reader.union('nhf')
            frame.loc[0, 'source_tier_highest'] = 2
            frame.to_csv(path, index=False)
            manifest_path = context.output / 'manifest.json'
            manifest = json.loads(manifest_path.read_text())
            manifest['outputs'][path.name]['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'tier disagrees'):
                reader.union('nhf')

if __name__ == '__main__':
    unittest.main()
