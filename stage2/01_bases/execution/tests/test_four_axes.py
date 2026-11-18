"""ADR0073 behavior through public cross-disease/axis/view boundaries."""
import json
import tempfile
import sys
import unittest
from pathlib import Path
import numpy as np

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / 'output/api'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bases import gene_axis


class FourAxisTests(unittest.TestCase):
    def test_f1000_uses_full_union_denominator_and_penalises_actual_deg_1001(self):
        ranks = np.full(1002, 1002, dtype=np.int64)
        ranks[:1000] = np.arange(1, 1001)
        result = gene_axis.cross_disease_rank_sum([ranks], [1001], categories=('nhf',))
        self.assertEqual(result['f1000'][0], 1)
        self.assertEqual(result['f1000'][1000], 0)
        self.assertAlmostEqual(result['rank_sum'][999], 1000 / 1001)
        self.assertAlmostEqual(result['rank_sum'][1000], 1002 / 1001)
        np.testing.assert_array_equal(gene_axis.f1000_count([ranks], [1001], categories=('nhf',)), result['f1000'])

    def test_four_exact_nested_axes_keep_hvg_half_and_record_fallback(self):
        self.assertEqual(gene_axis.AXIS_SIZES, (2000, 3000, 4000, 5000))
        hvg = np.arange(1, 6001)
        support = np.zeros(6000, dtype=np.int64)
        support[3500:4500] = 1
        support[4500:5000] = 2
        sums = np.ones(6000)
        previous = set()
        for size, expected_fallback in [(2000, 0), (3000, 0), (4000, 500), (5000, 1000)]:
            result = gene_axis.axis_fill(hvg, support, sums, size, n_diseases=2)
            genes = result['axis']
            self.assertEqual(len(genes), size)
            self.assertEqual(len(set(genes)), size)
            self.assertTrue(previous <= set(genes))
            self.assertTrue(set(range(size // 2)) <= set(genes))
            np.testing.assert_array_equal(genes, np.sort(genes))
            self.assertEqual(result['fallback'], expected_fallback)
            previous = set(genes)

    def test_axis_reader_rejects_mixed_algorithm_before_any_payload_read(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'manifest.json').write_text(json.dumps({'schema': gene_axis.SCHEMA, 'release': {'algorithm': 'f2000-deg2000-v2'}}))
            with self.assertRaisesRegex(ValueError, 'algorithm'):
                gene_axis.GeneAxis(root).manifest()

    def test_only_5k_expression_is_released_but_four_axes_remain(self):
        from bases import view
        self.assertEqual(tuple(view.VIEW_NAMES), (5000,))
        self.assertEqual(gene_axis.AXIS_SIZES, (2000, 3000, 4000, 5000))
        with tempfile.TemporaryDirectory() as temp:
            for size in (2000, 3000, 4000):
                with self.assertRaisesRegex(view.ViewError, 'released view size'):
                    view.build(Path(temp) / str(size), corpus_dir=Path(temp),
                               axis_dir=Path(temp), size=size)

    def test_view_reader_rejects_a_mixed_axis_algorithm(self):
        from bases import view
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'manifest.json').write_text(json.dumps({'schema': view.SCHEMA, 'axis': {'algorithm': 'retired-axis'}}))
            with self.assertRaisesRegex(view.ViewError, 'algorithm'):
                view.open_view(root).manifest()


if __name__ == '__main__':
    unittest.main()
