"""Tests for the cross-disease DEG fusion and the released axes (ADR 0071).

Everything is a hand-built 8-gene fixture except the real-file tests: the rule
document and the two committed figure PNGs are regenerated from their tracked
inputs and compared byte-for-byte, and the document's input wiring, its
source-table/per-disease contracts and its sample-size footnote are asserted
against the committed file. No release payload (h5ad, ``source/`` tree or
gitignored intermediate) is read.

The fusion, the ``f2000`` count, the axis fill rule (including the tenv-HVG
fallback), the readers and the composition/coverage products are exercised at
the public API / product boundary; one seam only (ADR 0071 / spec #200).
Determinism is tested by building the same fixture twice into scratch candidate
roots and comparing every digest. The retired ``pN`` has no column, field, reader
or test here: a regression that reintroduces one is meant to fail.
"""
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

EXECUTION = Path(__file__).resolve().parents[1]
STEP = EXECUTION.parent
sys.path.insert(0, str(EXECUTION))
sys.path.insert(0, str(STEP / 'output/api'))
sys.path.insert(0, str(STEP / 'figure/script'))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
import build_gene_axis as driver
import render_gene_axis_rank_sum_histogram as renderer
import render_gene_axis_source_topn_coverage as coverage_renderer
from bases import deg_rank, disease_deg, gene_axis
from bases.gene_axis import (
    # The fixtures exercise the *released* five-category products; the declared
    # roster is ``gene_axis.CATEGORIES`` (it gains ``af``, which #208 publishes
    # and moves into the release set).
    ALGORITHM, CANDIDATE_COLUMNS, CATEGORY_LABELS,
    CATEGORY_RANK_COLUMNS, COMPOSITION_COLUMNS, COVERAGE_COLUMNS, FILL_RULE,
    GeneAxis, RELEASE_CATEGORIES as CATEGORIES,
    axis_fill, candidate_columns, candidate_table, category_rank,
    cross_disease_rank_sum, deg2000_composition, f2000_count,
    f2000_distribution, source_topn_coverage,
)

AXIS_SYMBOLS = ['G0', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7']
AXIS_SIZE = len(AXIS_SYMBOLS)


def axis_fixture():
    return pd.DataFrame({
        'gene_index': np.arange(len(AXIS_SYMBOLS), dtype=np.int64),
        'symbol': AXIS_SYMBOLS,
        'hgnc_id': [f'HGNC:{index + 1}' for index in range(len(AXIS_SYMBOLS))],
        'ensembl_id': [f'ENSG0000000000{index + 1}' for index in range(len(AXIS_SYMBOLS))],
    })

class ReExportContractTests(unittest.TestCase):
    """#186 contract: gene_axis re-exports the DEG/rank functions, never copies them."""

    def test_deg_rank_contract_is_re_exported_by_identity(self):
        for name in ('SymbolResolver', 'complete_rank', 'deg_effect', 'per_gene_deg',
                     'rank_scores', 'ALIGNMENT_CATEGORIES', 'PADJ_FLOOR',
                     'PADJ_THRESHOLD', 'normalized_rank', 'AXIS_SIZE'):
            self.assertIs(getattr(gene_axis, name), getattr(deg_rank, name),
                          f'bases.gene_axis.{name} is not bases.deg_rank.{name}')



class LevelOneTests(unittest.TestCase):
    """Level 1 over the v2 DEG frames: sum(rank/n_DEG), ties by gene_index."""

    def test_single_dataset_category_is_its_own_order(self):
        ranks = np.array([8, 1, 7, 2, 6, 3, 5, 4], dtype=np.int64)
        result = category_rank([ranks], [8])
        np.testing.assert_array_equal(result['rank'], ranks)
        self.assertEqual(result['n_datasets'], 1)
        self.assertTrue(result['measured'].all())
        # ``deg`` is the v2 name for the same union mask; ``measured`` is the
        # retired ADR 0066 name a pre-v2 reader still imports.
        self.assertIs(result['deg'], result['measured'])

    def test_sum_of_normalized_ranks_and_gene_index_tie_break(self):
        # A = 1..8, B = 8,1,2,...,7 over 8 measured genes. Sum/8:
        # [9,3,5,7,9,11,13,15]/8 -> ascending 3,5,7,9,9,11,13,15 with the two
        # 9/8 ties broken by gene_index (gene 0 before gene 4).
        a = np.arange(1, 9, dtype=np.int64)
        b = np.array([8, 1, 2, 3, 4, 5, 6, 7], dtype=np.int64)
        result = category_rank([a, b], [8, 8])
        np.testing.assert_array_equal(result['rank'], np.array([4, 1, 2, 3, 5, 6, 7, 8]))
        np.testing.assert_allclose(
            result['sum'], np.array([9, 3, 5, 7, 9, 11, 13, 15]) / 8.0)

    def test_non_deg_contributes_n_deg_plus_one_over_n_deg(self):
        # Six DEGs rank 1..6; the two other genes carry n_DEG + 1 = 7, so their
        # rank_norm is 7/6 (ADR 0071 decision 3), never clamped to 1.0.
        ranks = np.array([1, 2, 3, 4, 5, 6, 7, 7], dtype=np.int64)
        result = category_rank([ranks], [6])
        np.testing.assert_allclose(result['sum'], [1/6, 2/6, 3/6, 4/6, 5/6, 1.0, 7/6, 7/6])
        np.testing.assert_array_equal(result['rank'], np.arange(1, 9))
        np.testing.assert_array_equal(
            result['deg'], [True, True, True, True, True, True, False, False])

    def test_zero_deg_source_is_rejected_as_non_informative(self):
        # ADR 0071 makes rank_norm undefined for n_DEG = 0; the caller must
        # exclude the source and record the exclusion (ADR 0068).
        with self.assertRaises(ValueError) as raised:
            category_rank([np.array([1, 1, 1, 1, 1, 1, 1, 1])], [0])
        self.assertIn('n_DEG = 0', str(raised.exception))

    def test_incomplete_deg_frame_fails_closed(self):
        # Two genes sit at n_DEG + 1 although n_DEG = 6 claims six DEGs.
        with self.assertRaises(ValueError):
            category_rank([np.array([1, 2, 3, 4, 5, 7, 7, 7])], [6])
        # A rank beyond n_DEG + 1.
        with self.assertRaises(ValueError):
            category_rank([np.array([1, 2, 3, 4, 5, 6, 7, 9])], [6])

    def test_level_one_rank_is_a_complete_permutation(self):
        a = np.arange(1, 9, dtype=np.int64)
        b = np.array([2, 1, 4, 3, 6, 5, 8, 7], dtype=np.int64)
        result = category_rank([a, b], [8, 8])
        self.assertEqual(sorted(result['rank'].tolist()), list(range(1, 9)))
        # every pair ties at 3,7,11,15 / 8, so gene_index orders each pair
        np.testing.assert_array_equal(result['rank'], np.arange(1, 9))

    def test_out_of_contract_rank_fails_closed(self):
        with self.assertRaises(ValueError):
            category_rank([np.array([1, 2, 3, 4, 5, 6, 7, 9])], [8])
        with self.assertRaises(ValueError):
            category_rank([np.array([1, 2, 3, 4, 5, 6, 7, 8])], [0])
        with self.assertRaises(ValueError):
            category_rank([], [])



class CrossDiseaseFusionTests(unittest.TestCase):
    """Level 2: ``f2000`` counts diseases; ``rank_sum`` divides by the union size."""

    @staticmethod
    def frame(n_deg, ranks):
        frame = np.full(AXIS_SIZE, n_deg + 1, dtype=np.int64)
        for gene, rank in ranks.items():
            frame[gene] = rank
        return frame

    def test_f2000_counts_diseases_not_datasets(self):
        a = self.frame(2, {0: 1, 1: 2})
        b = self.frame(2, {1: 1, 2: 2})
        c = self.frame(1, {2: 1})
        result = cross_disease_rank_sum([a, b, c], [2, 2, 1],
                                        categories=('nhf', 'icm', 'af'))
        np.testing.assert_array_equal(result['f2000'], [1, 2, 2, 0, 0, 0, 0, 0])
        np.testing.assert_array_equal(
            f2000_count([a, b, c], [2, 2, 1], categories=('nhf', 'icm', 'af')),
            result['f2000'])
        self.assertEqual(int(result['f2000'].max()), 2)
        # one mask per active disease, in the requested order
        self.assertEqual(len(result['deg2000']), 3)
        np.testing.assert_array_equal(result['deg2000'][0], [True, True] + [False] * 6)

    def test_rank_sum_uses_the_union_denominator_and_penalises_outside_genes(self):
        a = self.frame(3, {0: 1, 1: 2, 2: 3})   # union 3, every member ranked
        b = self.frame(2, {0: 2, 1: 1})         # union 2, every member ranked
        result = cross_disease_rank_sum([a, b], [3, 2], categories=('nhf', 'icm'))
        outside = 4 / 3 + 3 / 2                    # (n_d + 1) / n_d for both diseases
        np.testing.assert_allclose(
            result['rank_sum'],
            [1 / 3 + 2 / 2, 2 / 3 + 1 / 2, 3 / 3 + 3 / 2,
             outside, outside, outside, outside, outside])
        # A gene outside DEG2000_d contributes (n_d + 1) / n_d, never 0 and never
        # an omitted term: every disease adds its own penalty.
        self.assertGreater(result['rank_sum'][3], 4 / 3)

    def test_rank_sum_denominator_is_the_union_size_not_the_deg2000_cut(self):
        # The published frame is the DEG2000 window; the denominator is the
        # disease's axis-resolved DEG-union size, so a union larger than the cut
        # is expressible with a smaller limit.
        frame = self.frame(3, {0: 1, 1: 2})
        result = cross_disease_rank_sum([frame], [3], categories=('nhf',), limit=2)
        np.testing.assert_array_equal(result['f2000'], [1, 1, 0, 0, 0, 0, 0, 0])
        np.testing.assert_allclose(result['rank_sum'][:3],
                                   [1 / 3, 2 / 3, 4 / 3])
        # A member of the union that the 2,000-cap cut out must sit at n_deg + 1,
        # not at its within-union rank.
        bad = np.full(AXIS_SIZE, 4, dtype=np.int64)
        bad[[0, 1, 2]] = [1, 2, 3]
        with self.assertRaisesRegex(ValueError, 'n_deg\\+1'):
            cross_disease_rank_sum([bad], [3], categories=('nhf',), limit=2)

    def test_tenv_is_not_a_rank_sum_term(self):
        frame = self.frame(3, {0: 1, 1: 2, 2: 3})
        axis = axis_fixture()
        forward = np.arange(1, AXIS_SIZE + 1, dtype=np.int64)
        reverse = np.arange(AXIS_SIZE, 0, -1, dtype=np.int64)
        one = candidate_table(axis, forward, [frame], [3], categories=('nhf',))
        two = candidate_table(axis, reverse, [frame], [3], categories=('nhf',))
        np.testing.assert_allclose(one.rank_sum, two.rank_sum)
        np.testing.assert_array_equal(one.f2000, two.f2000)
        self.assertFalse(np.array_equal(one.tenv_hvg_rank, two.tenv_hvg_rank))

    def test_out_of_contract_frames_fail_closed(self):
        good = self.frame(3, {0: 1, 1: 2, 2: 3})
        with self.assertRaises(ValueError):
            cross_disease_rank_sum([], [])                       # no disease
        with self.assertRaises(ValueError):
            cross_disease_rank_sum([good], [])                   # missing n_deg
        with self.assertRaises(ValueError):
            cross_disease_rank_sum([good], [3], categories=('nhf', 'icm'))
        with self.assertRaises(ValueError):
            cross_disease_rank_sum([good, good], [3, 3],
                                   categories=('nhf', 'nhf'))    # repeated
        with self.assertRaises(ValueError):
            cross_disease_rank_sum([good], [3], categories=('unknown',))
        with self.assertRaises(ValueError):
            cross_disease_rank_sum([good[:4]], [3])              # wrong shape
        with self.assertRaises(ValueError):
            cross_disease_rank_sum([self.frame(3, {0: 1, 1: 2})], [0])
        # A frame whose DEG2000 ranks are not a complete 1..k assignment.
        with self.assertRaises(ValueError):
            cross_disease_rank_sum([self.frame(3, {0: 2, 1: 3})], [3])

    def test_fusion_is_deterministic(self):
        frames = [self.frame(3, {0: 1, 1: 2, 2: 3}),
                  self.frame(2, {2: 1, 3: 2})]
        one = cross_disease_rank_sum(frames, [3, 2], categories=('nhf', 'icm'))
        two = cross_disease_rank_sum(frames, [3, 2], categories=('nhf', 'icm'))
        np.testing.assert_array_equal(one['rank_sum'], two['rank_sum'])
        np.testing.assert_array_equal(one['f2000'], two['f2000'])


class CandidateTableTests(unittest.TestCase):
    """The master table carries ``f2000`` and no retired ``pN`` column."""

    @staticmethod
    def frame(n_deg, ranks):
        return CrossDiseaseFusionTests.frame(n_deg, ranks)

    def test_columns_rows_and_values(self):
        axis = axis_fixture()
        tenv = np.array([8, 7, 6, 5, 4, 3, 2, 1], dtype=np.int64)
        nhf = self.frame(3, {0: 1, 1: 2, 2: 3})
        icm = self.frame(2, {4: 1, 0: 2})
        table = candidate_table(axis, tenv, [nhf, icm], [3, 2],
                                categories=('nhf', 'icm'))
        self.assertEqual(list(table.columns), list(candidate_columns()))
        self.assertNotIn('p2000', table.columns)
        self.assertNotIn('p5000', table.columns)
        self.assertNotIn('p2000', list(CANDIDATE_COLUMNS))
        self.assertEqual(len(table), AXIS_SIZE)
        self.assertEqual(table.gene_index.tolist(), list(range(AXIS_SIZE)))
        self.assertEqual(table.tenv_hvg_rank.tolist(), tenv.tolist())
        # DEG2000 ranks are published, non-members are empty.
        self.assertEqual(table.nhf_rank.iloc[:3].tolist(), [1.0, 2.0, 3.0])
        self.assertTrue(table.nhf_rank.iloc[3:].isna().all())
        self.assertEqual(table.icm_rank.iloc[4], 1.0)
        self.assertEqual(table.icm_rank.iloc[0], 2.0)
        self.assertTrue(table.icm_rank.iloc[1:4].isna().all())
        np.testing.assert_array_equal(table.f2000, [2, 1, 1, 0, 1, 0, 0, 0])
        np.testing.assert_allclose(table.rank_sum[0], 1 / 3 + 2 / 2)
        np.testing.assert_allclose(table.rank_sum[1], 2 / 3 + 3 / 2)

    def test_excluded_category_keeps_an_empty_column(self):
        axis = axis_fixture()
        tenv = np.arange(1, AXIS_SIZE + 1, dtype=np.int64)
        table = candidate_table(axis, tenv, [self.frame(2, {0: 1, 1: 2})], [2],
                                categories=('nhf',))
        for category in CATEGORIES:
            self.assertIn(f'{category}_rank', table.columns)
        self.assertTrue(table.icm_rank.isna().all())
        self.assertEqual(int(table.f2000.max()), 1)

    def test_axis_order_and_tenv_permutation_are_required(self):
        axis = axis_fixture().iloc[::-1].reset_index(drop=True)
        frame = self.frame(2, {0: 1, 1: 2})
        with self.assertRaises(ValueError):
            candidate_table(axis, np.arange(1, 9), [frame], [2],
                            categories=('nhf',))
        with self.assertRaises(ValueError):
            candidate_table(axis_fixture(),
                            np.array([1, 1, 2, 3, 4, 5, 6, 7]), [frame], [2],
                            categories=('nhf',))


class AxisFillTests(unittest.TestCase):
    """The declared ``0.5*N + 0.5*N`` fill, its ordering and its fallback."""

    TENV = np.array([8, 7, 6, 5, 4, 3, 2, 1], dtype=np.int64)

    def fill(self, f2000, rank_sum, n, n_diseases=3):
        return axis_fill(self.TENV, np.asarray(f2000, dtype=np.int64),
                         np.asarray(rank_sum, dtype=np.float64), n=n,
                         n_diseases=n_diseases)

    def test_fill_takes_exactly_half_tenv_and_half_the_disease_order(self):
        f2000 = [0, 2, 3, 1, 2, 1, 0, 0]
        rank_sum = [7.0, 5.0, 4.0, 6.0, 5.5, 6.5, 7.5, 7.5]
        result = self.fill(f2000, rank_sum, 4)
        np.testing.assert_array_equal(result['axis'], [1, 2, 6, 7])
        self.assertEqual(result['tenv_half'], 2)
        self.assertEqual(result['from_deg2000'], 2)
        self.assertEqual(result['fallback'], 0)
        # G_6 keeps the same order: the disease half is a prefix, so the axes are
        # nested when the fallback is unused.
        larger = self.fill(f2000, rank_sum, 6)
        np.testing.assert_array_equal(larger['axis'], [1, 2, 4, 5, 6, 7])
        self.assertEqual(larger['from_deg2000'], 3)
        self.assertEqual(larger['fallback'], 0)

    def test_f2000_ties_break_on_rank_sum(self):
        f2000 = [0, 0, 1, 0, 1, 1, 0, 0]
        # All three disease candidates share f2000 = 1; the two best rank_sums win.
        result = self.fill(f2000, [9, 9, 2.0, 9, 1.0, 1.5, 9, 9], 4)
        np.testing.assert_array_equal(result['axis'], [4, 5, 6, 7])
        self.assertEqual(result['from_deg2000'], 2)
        self.assertNotIn(2, result['axis'].tolist())

    def test_rank_sum_ties_break_on_gene_index(self):
        f2000 = [0, 0, 1, 0, 1, 1, 0, 0]
        # Genes 2 and 4 tie on both keys at the boundary; gene_index orders them.
        at_boundary = self.fill(f2000, [9, 9, 1.0, 9, 2.0, 1.0, 9, 9], 4)
        np.testing.assert_array_equal(at_boundary['axis'], [2, 5, 6, 7])
        asymmetric = self.fill(f2000, [9, 9, 1.0, 9, 1.0, 2.0, 9, 9], 4)
        np.testing.assert_array_equal(asymmetric['axis'], [2, 4, 6, 7])

    def test_exact_half_boundary_needs_no_fallback(self):
        f2000 = [0, 0, 1, 0, 0, 1, 0, 0]
        result = self.fill(f2000, [9, 9, 3.0, 9, 9, 4.0, 9, 9], 4)
        self.assertEqual(result['from_deg2000'], 2)
        self.assertEqual(result['fallback'], 0)
        np.testing.assert_array_equal(result['axis'], [2, 5, 6, 7])

    def test_fallback_uses_unused_tenv_genes_in_tenv_rank_order(self):
        f2000 = [0, 0, 1, 0, 0, 0, 0, 0]
        result = self.fill(f2000, [9, 9, 1.0, 9, 9, 9, 9, 9], 4)
        # The only disease candidate is gene 2; the missing slot comes from the
        # best-ranked unused tenv gene, gene 5 (tenv rank 3).
        self.assertEqual(result['from_deg2000'], 1)
        self.assertEqual(result['fallback'], 1)
        np.testing.assert_array_equal(result['axis'], [2, 5, 6, 7])
        larger = self.fill(f2000, [9, 9, 1.0, 9, 9, 9, 9, 9], 6)
        self.assertEqual(larger['from_deg2000'], 1)
        self.assertEqual(larger['fallback'], 2)
        np.testing.assert_array_equal(larger['axis'], [2, 3, 4, 5, 6, 7])

    def test_empty_deg2000_union_fills_entirely_from_tenv(self):
        result = self.fill([0] * AXIS_SIZE, [9.0] * AXIS_SIZE, 4, n_diseases=0)
        self.assertEqual(result['from_deg2000'], 0)
        self.assertEqual(result['fallback'], 2)
        np.testing.assert_array_equal(result['axis'], [4, 5, 6, 7])

    def test_fallback_can_exhaust_every_unused_tenv_gene(self):
        # N = the frozen axis and no DEG2000 at all: the second half consumes every
        # unused tenv-HVG gene, so the axis is the whole frozen space.
        result = self.fill([0] * AXIS_SIZE, [9.0] * AXIS_SIZE, AXIS_SIZE,
                           n_diseases=0)
        self.assertEqual(result['from_deg2000'], 0)
        self.assertEqual(result['fallback'], AXIS_SIZE // 2)
        np.testing.assert_array_equal(result['axis'], np.arange(AXIS_SIZE))

    def test_a_gene_without_disease_evidence_never_enters_the_disease_half(self):
        f2000 = [0, 0, 1, 0, 0, 0, 0, 0]
        rank_sum = [0.1, 9, 9.9, 9, 9, 9, 9, 9]
        # Gene 0 has the best rank_sum but no disease evidence: it is not a
        # disease-half candidate, and the disease half's fallback is tenv-ordered
        # (gene 5), so gene 0 is excluded even though its tenv rank is last.
        result = self.fill(f2000, rank_sum, 4)
        self.assertEqual(result['from_deg2000'], 1)
        self.assertEqual(result['fallback'], 1)
        self.assertNotIn(0, result['axis'].tolist())

    def test_fill_is_deterministic(self):
        f2000 = [0, 1, 1, 2, 0, 1, 0, 0]
        rank_sum = [5, 4, 4, 1, 5, 4, 5, 5]
        one = self.fill(f2000, rank_sum, 6)
        two = self.fill(f2000, rank_sum, 6)
        np.testing.assert_array_equal(one['axis'], two['axis'])
        self.assertEqual(one['from_deg2000'], two['from_deg2000'])
        self.assertEqual(one['fallback'], two['fallback'])

    def test_out_of_contract_fill_fails_closed(self):
        f2000 = [0, 1, 1, 2, 0, 1, 0, 0]
        rank_sum = [5, 4, 4, 1, 5, 4, 5, 5]
        with self.assertRaises(ValueError):
            self.fill(f2000, rank_sum, 7)                 # odd N
        with self.assertRaises(ValueError):
            self.fill(f2000, rank_sum, 10)                # larger than the axis
        with self.assertRaises(ValueError):
            axis_fill(self.TENV[:7], f2000[:7], rank_sum[:7], n=4)
        with self.assertRaises(ValueError):
            axis_fill(np.array([1, 1, 2, 3, 4, 5, 6, 7]), f2000, rank_sum, n=4)
        with self.assertRaises(ValueError):
            self.fill([0, 3, 0, 0, 0, 0, 0, 0], rank_sum, 4,
                      n_diseases=2)   # f2000 > diseases
        with self.assertRaises(ValueError):
            self.fill([0, -1, 0, 0, 0, 0, 0, 0], rank_sum, 4)
        with self.assertRaises(ValueError):
            self.fill(f2000, [np.nan] * AXIS_SIZE, 4)


class HistogramAndInflectionTests(unittest.TestCase):
    """Evidence panel: rank-sum histogram and the diagnostic-only inflection."""

    @staticmethod
    def bent_curve(n=1000, knee=300):
        # piecewise linear: y = x up to the knee, then a tenth of that slope.
        # Both criteria then point at the same position, as the design intends.
        x = np.arange(1, n + 1, dtype=np.float64)
        y = np.where(x <= knee, x, knee + 0.1 * (x - knee))
        return x, y

    def test_histogram_covers_every_candidate_and_keeps_empty_bins(self):
        values = np.array([0.1, 0.2, 0.2, 6.9])
        frame = gene_axis.rank_sum_histogram(values, low=0.0, high=7.0, bins=7)
        self.assertEqual(int(frame['count'].sum()), len(values))
        self.assertEqual(frame['count'].tolist(), [3, 0, 0, 0, 0, 0, 1])
        np.testing.assert_allclose(frame['bin_left'], np.arange(7, dtype=float))
        np.testing.assert_allclose(frame['bin_center'], np.arange(7) + 0.5)

    def test_histogram_rejects_values_outside_the_declared_interval(self):
        with self.assertRaises(ValueError):
            gene_axis.rank_sum_histogram(np.array([-0.1, 1.0]), low=0.0, high=7.0)
        with self.assertRaises(ValueError):
            gene_axis.rank_sum_histogram(np.array([7.5]), low=0.0, high=7.0)

    def test_synthetic_knee_is_detected(self):
        x, y = self.bent_curve()
        result = gene_axis.detect_inflection(x, y, min_n=10)
        self.assertTrue(result['inflection'])
        self.assertIsNotNone(result['n'])
        self.assertGreater(result['n'], 250)
        self.assertLess(result['n'], 350)
        self.assertIn('parameters', result)
        self.assertIn('sensitivity', result['parameters'])
        self.assertIsNotNone(result['kneedle'])
        self.assertIsNotNone(result['second_difference'])

    def test_straight_line_has_no_inflection(self):
        x = np.arange(1, 1001, dtype=np.float64)
        y = 2.0 + 0.001 * x
        result = gene_axis.detect_inflection(x, y, min_n=10)
        self.assertFalse(result['inflection'])
        self.assertIsNone(result['n'])
        self.assertIsNotNone(result['parameters'])

    def test_stability_reports_leave_one_source_out_agreement(self):
        x, y = self.bent_curve()
        stable = gene_axis.inflection_stability([(x, y), (x, y * 1.001), (x, y * 0.999)])
        self.assertTrue(stable['stable'])
        self.assertEqual(stable['n_replays'], 2)
        self.assertEqual(stable['n_same_decision'], 2)
        straight = 2.0 + 0.001 * x
        unstable = gene_axis.inflection_stability([(x, y), (x, straight)])
        self.assertFalse(unstable['stable'])
        self.assertEqual(unstable['n_same_decision'], 0)





# ---- hermetic release fixture ------------------------------------------------

#: The declared released axis sizes of the fixture (halves 2 and 3 genes).
FIXTURE_AXIS_SIZES = (4, 6)
#: tenv-HVG ranks: gene 7 is the strongest (rank 1), gene 0 the weakest.
FIXTURE_TENV_RANK = np.array([8, 7, 6, 5, 4, 3, 2, 1], dtype=np.int64)

#: The default fixture: the six released diseases, each a single informative
#: source whose published ``DEG2000_d`` is ``(n_deg, {gene_index: disease_rank})``.
#: ``af`` covers the tenv-first-half gene 7 only, so it adds the constant
#: ``(n_af + 1) / n_af = 2`` to every non-AF gene's rank_sum and leaves the
#: second-half order unchanged. The level-2 vectors are hand-computable:
#: f2000 = [2, 2, 3, 2, 2, 1, 0, 1] and
#: rank_sum = [8, 47/6, 7.5, 8.5, 49/6, 8.5, 9.5, 8.5] (gene 7 gains AF's 1/1
#: because it is AF's only DEG2000 gene, and still pays every other disease's
#: (n_d + 1) / n_d outside penalty).
FIXTURE_DEG2000 = {
    'nhf': (3, {0: 1, 1: 2, 2: 3}),
    'icm': (3, {1: 1, 2: 2, 3: 3}),
    'af': (1, {7: 1}),
    'genetic_cardiomyopathy': (3, {2: 1, 3: 2, 4: 3}),
    'cardiotox_TOP2': (2, {4: 1, 0: 2}),
    'cardiotox_other': (1, {5: 1}),
}


def published_frame(n_deg, ranks):
    """A complete frozen-axis DEG2000 frame: rank inside, ``n_deg + 1`` outside."""
    frame = np.full(AXIS_SIZE, n_deg + 1, dtype=np.int64)
    for gene, rank in ranks.items():
        frame[gene] = rank
    return frame


def fusion_disease_record(directory, n_deg, ranks):
    """One valid published disease record of the ``disease-fusion`` catalog."""
    label = CATEGORY_LABELS[directory]
    rows = pd.DataFrame({
        'gene_index': np.array(list(ranks), dtype=np.int64),
        'disease_rank': np.array(list(ranks.values()), dtype=np.int64),
    })
    rows = rows.loc[:, list(disease_deg.DEG2000_COLUMNS[:2])]
    rows['disease_rank_sum'] = rows.disease_rank / 10.0
    rows = rows.loc[:, list(disease_deg.DEG2000_COLUMNS)]
    deg2000 = min(int(n_deg), int(disease_deg.DEG2000_LIMIT)) if n_deg else 0
    record = {
        'category': label, 'directory': directory, 'published': True, 'tier': 1,
        'evidence_type': 'cm_pseudobulk',
        'declared_source_count': 0 if n_deg == 0 else 1,
        'informative_source_count': 0 if n_deg == 0 else 1,
        'n_deg': int(n_deg), 'deg2000': int(deg2000),
        'deg2000_file': f'{directory}.deg2000.csv',
        'sources': ([] if n_deg == 0 else [{
            'dataset_id': f'{directory}_a', 'tier': 1, 'informative': True,
            'n_deg': int(n_deg), 'n_deg_in_deg2000': int(deg2000)}]),
        'excluded_sources': [], 'not_published': [],
    }
    return record, rows


def fusion_catalog(genes_path, records):
    """The published ``disease-fusion`` catalog over the fixture axis."""
    return {
        'schema': disease_deg.FUSION_SCHEMA,
        'axis': {'path': 'stage2/01_bases/output/gene-identity/genes.csv',
                 'sha256': hashlib.sha256(Path(genes_path).read_bytes()).hexdigest(),
                 'size': AXIS_SIZE},
        'limit': int(disease_deg.DEG2000_LIMIT),
        'diseases': records,
    }


def write_fixture(root, *, deg2000=None):
    """Write a hermetic mini repository: the frozen axis, tenv-HVG and the fusion.

    ``deg2000`` maps a released category to ``(n_deg, {gene_index: rank})``; a
    category omitted from it is published with an empty ``DEG2000`` (the
    representable zero-union disease). The default is :data:`FIXTURE_DEG2000`.
    """
    root = Path(root)
    step = root / 'stage2/01_bases'
    genes = axis_fixture()
    genes_path = step / 'output/gene-identity/genes.csv'
    genes_path.parent.mkdir(parents=True, exist_ok=True)
    genes.to_csv(genes_path, index=False, lineterminator='\n')
    tenv = pd.DataFrame({
        'gene_index': np.arange(AXIS_SIZE),
        'log_mean': np.zeros(AXIS_SIZE),
        'log_dispersion': np.zeros(AXIS_SIZE),
        'normalized_dispersion': -FIXTURE_TENV_RANK.astype(float),
        'rank': FIXTURE_TENV_RANK,
    })
    tenv_path = step / 'output/tenv-hvg/descriptive.csv'
    tenv_path.parent.mkdir(parents=True, exist_ok=True)
    tenv.to_csv(tenv_path, index=False, lineterminator='\n')
    fusion_root = step / 'output/disease-fusion'
    fusion_root.mkdir(parents=True, exist_ok=True)
    definitions = dict(FIXTURE_DEG2000 if deg2000 is None else deg2000)
    records, outputs = [], {}
    for category in CATEGORIES:
        n_deg, ranks = definitions.get(category, (0, {}))
        record, rows = fusion_disease_record(category, n_deg, ranks)
        path = fusion_root / record['deg2000_file']
        rows.to_csv(path, index=False, lineterminator='\n')
        outputs[record['deg2000_file']] = {
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
        records.append(record)
    disease_deg.write_catalog(fusion_root / 'catalog.json',
                              fusion_catalog(genes_path, records))
    (fusion_root / 'manifest.json').write_text(json.dumps({
        'schema': disease_deg.FUSION_RELEASE_SCHEMA,
        'outputs': dict(sorted(outputs.items())),
    }, indent=2, sort_keys=True) + '\n')
    return driver.default_config(
        root=root,
        output=step / 'output/gene-axis',
        figure_data=step / 'figure/data/gene_axis_rank_sum_histogram.csv',
        axis_sizes=FIXTURE_AXIS_SIZES)


def rewrite_fusion_catalog(config, mutate):
    """Rewrite the fixture fusion catalog through a mutation (re-signing it)."""
    path = Path(config.fusion_root) / 'catalog.json'
    catalog = json.loads(path.read_text())
    catalog.pop('catalog_sha256', None)
    mutate(catalog)
    disease_deg.write_catalog(path, catalog)
    return catalog


def rewrite_fusion_output(config, name, mutate):
    """Rewrite one fixture DEG2000 payload and re-pin it in the fusion manifest."""
    path = Path(config.fusion_root) / name
    frame = pd.read_csv(path)
    frame = mutate(frame)
    frame.to_csv(path, index=False, lineterminator='\n')
    manifest_path = Path(config.fusion_root) / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['outputs'][name] = {
        'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')

def manifest_strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from manifest_strings(key)
            yield from manifest_strings(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from manifest_strings(item)



class DeclaredRosterTests(unittest.TestCase):
    """The six-disease roster is the released roster after #208.

    ``gene_axis.UNPUBLISHED_CATEGORIES`` is empty, so ``RELEASE_CATEGORIES`` is the
    full six-disease ``CATEGORIES`` roster and the driver fails closed on a
    missing released product or a published disease outside the roster.
    """

    def test_all_six_declared_categories_are_released(self):
        self.assertEqual(gene_axis.UNPUBLISHED_CATEGORIES, ())
        self.assertEqual(gene_axis.RELEASE_CATEGORIES, gene_axis.CATEGORIES)
        self.assertIn('af', gene_axis.RELEASE_CATEGORIES)

    def test_the_release_contract_carries_the_six_released_categories(self):
        self.assertEqual(len(gene_axis.CATEGORY_RANK_COLUMNS), 6)
        self.assertIn('af_rank', gene_axis.CATEGORY_RANK_COLUMNS)
        # The retired pN columns have no successor in the v2 contract.
        self.assertNotIn('p2000', list(gene_axis.CANDIDATE_COLUMNS))
        self.assertNotIn('p5000', list(gene_axis.CANDIDATE_COLUMNS))
        self.assertFalse(hasattr(gene_axis, 'fusion_rank_sum'))
        self.assertFalse(hasattr(gene_axis, 'pn_distribution'))
        self.assertFalse(hasattr(gene_axis, 'winning_tier'))

    def test_the_six_category_release_builds(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            result = driver.publish(config)
            self.assertEqual(result['result'], 'PASS')
            self.assertEqual(result['n_diseases'], 6)
            self.assertEqual(driver.verify(config)['result'], 'PASS')

    def test_a_missing_released_disease_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            rewrite_fusion_catalog(
                config, lambda catalog: catalog.__setitem__(
                    'diseases', [record for record in catalog['diseases']
                                 if record['directory'] != 'icm']))
            with self.assertRaisesRegex(ValueError, 'missing'):
                driver.publish(config)

    def test_an_undeclared_fusion_directory_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))

            def add_unknown(catalog):
                extra = json.loads(json.dumps(catalog['diseases'][0]))
                extra['directory'] = 'unknown_category'
                extra['category'] = 'unknown_category'
                extra['deg2000_file'] = 'unknown_category.deg2000.csv'
                catalog['diseases'].append(extra)

            rewrite_fusion_catalog(config, add_unknown)
            with self.assertRaisesRegex(ValueError, 'undeclared'):
                driver.publish(config)

    def test_a_different_frozen_axis_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            rewrite_fusion_catalog(
                config, lambda catalog: catalog['axis'].__setitem__('sha256', '0' * 64))
            with self.assertRaisesRegex(ValueError, 'different frozen-axis digest'):
                driver.publish(config)


class ReleaseBuildTests(unittest.TestCase):
    """Candidate-then-replace publication, replay determinism and the reader seam."""

    def test_publish_products_and_reader_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            result = driver.publish(config)
            self.assertEqual(result['algorithm'], ALGORITHM)
            self.assertEqual(result['n_candidates'], AXIS_SIZE)
            self.assertEqual(result['fallback_second_half_slots'], {'4': 0, '6': 0})
            reader = GeneAxis(config.output)
            manifest = reader.manifest()
            self.assertEqual(manifest['schema'], gene_axis.SCHEMA)
            self.assertEqual(manifest['release']['algorithm'], ALGORITHM)
            self.assertEqual(manifest['release']['fill_rule'], FILL_RULE)
            self.assertIn('f2000', manifest['release'])
            self.assertIn('rank_sum', manifest['release'])
            self.assertEqual(reader.axis_sizes(), FIXTURE_AXIS_SIZES)
            table = reader.candidate_table()
            self.assertEqual(list(table.columns),
                             list(gene_axis.candidate_columns()))
            self.assertEqual(len(table), AXIS_SIZE)
            self.assertNotIn('p2000', table.columns)
            self.assertNotIn('p5000', table.columns)
            np.testing.assert_array_equal(table.f2000, [2, 2, 3, 2, 2, 1, 0, 1])
            np.testing.assert_allclose(
                table.rank_sum,
                [8.0, 47 / 6, 7.5, 8.5, 49 / 6, 8.5, 9.5, 8.5])
            np.testing.assert_array_equal(reader.axis(4), [1, 2, 6, 7])
            np.testing.assert_array_equal(reader.axis(6), [0, 1, 2, 5, 6, 7])
            self.assertEqual(manifest['axes']['4']['second_half'],
                             {'from_deg2000': 2, 'tenv_hvg_fallback': 0,
                              'n_deg2000_candidates': 6})
            self.assertEqual(manifest['axes']['6']['second_half']['from_deg2000'], 3)
            self.assertEqual(manifest['intersection']['n_intersection'], 4)
            self.assertEqual(manifest['intersection']['sizes'], [4, 6])
            self.assertEqual(manifest['fusion']['n_diseases'], 6)
            self.assertEqual(manifest['fusion']['diseases'], list(CATEGORIES))
            self.assertEqual(manifest['fusion']['excluded_diseases'], [])
            self.assertEqual(manifest['categories']['nhf']['n_deg'], 3)
            self.assertEqual(manifest['categories']['nhf']['deg2000'], 3)
            self.assertEqual(manifest['fusion']['n_declared_diseases'], 6)
            self.assertEqual(manifest['counts']['f2000_distribution'],
                             {'0': 1, '1': 2, '2': 4, '3': 1, '4': 0, '5': 0,
                              '6': 0})
            # The composition carries the declared halves, the DEG2000 numerator
            # and denominator, and the f2000 distribution; no pN record.
            composition = reader.composition()
            self.assertEqual(list(composition.columns), list(COMPOSITION_COLUMNS))
            self.assertEqual(len(composition), 2 * (3 + len(CATEGORIES)
                                                   + len(CATEGORIES) + 1))
            self.assertEqual(set(composition.record), {'axis_fill', 'deg2000', 'f2000'})
            axis2000 = composition[composition.axis.astype(str) == '4']
            fallback = axis2000[(axis2000.record == 'axis_fill')
                                & (axis2000.group == 'tenv_hvg_fallback')]
            self.assertEqual(int(fallback.n_genes.iloc[0]), 0)
            counts = axis2000[axis2000.record == 'f2000'].set_index('group')['n_genes']
            self.assertEqual(int(counts.sum()), 4)
            self.assertEqual(int(counts['0']), 1)   # gene 6 (tenv half only)
            self.assertEqual(int(counts['1']), 1)   # gene 7 (tenv half only)
            self.assertEqual(int(counts['2']), 1)   # gene 1
            self.assertEqual(int(counts['3']), 1)   # gene 2
            # Coverage: explicit numerator/denominator per axis and group.
            coverage = reader.coverage()
            self.assertEqual(list(coverage.columns), list(COVERAGE_COLUMNS))
            self.assertEqual(len(coverage), 2 * (1 + len(CATEGORIES)))
            tenv = coverage[coverage.group == 'tenv_hvg']
            self.assertEqual(tenv.denominator.tolist(), [4, 6])
            self.assertEqual(tenv.numerator.tolist(), [2, 4])
            self.assertEqual(tenv.rate.tolist(), [50.0, 66.67])
            nhf = coverage[coverage.group == 'nhf']
            self.assertEqual(nhf.denominator.tolist(), [3, 3])
            self.assertEqual(nhf.numerator.tolist(), [2, 3])
            # Files and absolute paths.
            self.assertIn('source_topn_coverage.csv',
                          ' '.join(manifest['outputs'].keys()))
            self.assertIn('gene_axis_source_topn_coverage.csv',
                          ' '.join(manifest['outputs'].keys()))
            coverage_figure = pd.read_csv(config.source_coverage_figure)
            self.assertEqual(list(coverage_figure.columns),
                             list(driver.SOURCE_COVERAGE_FIGURE_COLUMNS))
            self.assertEqual(len(coverage_figure),
                             2 + len(FIXTURE_AXIS_SIZES) * (len(CATEGORIES) + 1))
            self.assertEqual(len(coverage),
                             2 * (1 + len(CATEGORIES)))
            bars = coverage_figure[coverage_figure.record == 'bar']
            self.assertTrue((bars['count'] % 1 == 0).all())
            self.assertTrue(bars['rate'].between(0.0, 100.0).all())
            for text in manifest_strings(manifest):
                self.assertFalse(text.startswith('/'), f'absolute path in manifest: {text}')
            serialised = json.dumps(manifest)
            self.assertNotIn('p2000', serialised)
            self.assertNotIn('p5000', serialised)
            self.assertNotIn('"pN"', serialised)
            self.assertNotIn('n_sources_measured', serialised)

    def test_second_half_fallback_is_recorded_in_the_manifest(self):
        empty = {category: (0, {}) for category in CATEGORIES}
        empty['nhf'] = (1, {2: 1})
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary), deg2000=empty)
            result = driver.publish(config)
            self.assertEqual(result['fallback_second_half_slots'], {'4': 1, '6': 2})
            manifest = GeneAxis(config.output).manifest()
            self.assertEqual(manifest['axes']['4']['second_half'],
                             {'from_deg2000': 1, 'tenv_hvg_fallback': 1,
                              'n_deg2000_candidates': 1})
            self.assertEqual(manifest['axes']['6']['second_half'],
                             {'from_deg2000': 1, 'tenv_hvg_fallback': 2,
                              'n_deg2000_candidates': 1})
            reader = GeneAxis(config.output)
            np.testing.assert_array_equal(reader.axis(4), [2, 5, 6, 7])
            np.testing.assert_array_equal(reader.axis(6), [2, 3, 4, 5, 6, 7])
            composition = reader.composition()
            for size, fallback in (('4', 1), ('6', 2)):
                rows = composition[(composition.axis.astype(str) == size)
                                   & (composition.record == 'axis_fill')
                                   & (composition.group == 'tenv_hvg_fallback')]
                self.assertEqual(int(rows.n_genes.iloc[0]), fallback)
            # A fallback gene carries f2000 = 0; no disease gene was displaced.
            table = reader.candidate_table()
            self.assertEqual(int(table.f2000[2]), 1)
            self.assertEqual(int((table.f2000 == 0).sum()), 7)

    def test_excluded_disease_is_dropped_from_the_fusion_and_the_products(self):
        empty = {category: (0, {}) for category in CATEGORIES}
        empty['nhf'] = (2, {0: 1, 1: 2})
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary), deg2000=empty)
            driver.publish(config)
            reader = GeneAxis(config.output)
            manifest = reader.manifest()
            self.assertEqual(manifest['fusion']['n_diseases'], 1)
            self.assertEqual(manifest['fusion']['excluded_diseases'],
                             ['icm', 'af', 'genetic_cardiomyopathy',
                              'cardiotox_TOP2', 'cardiotox_other'])
            self.assertTrue(manifest['categories']['icm']['excluded_from_fusion'])
            self.assertIn('excluded_reason', manifest['categories']['icm'])
            table = reader.candidate_table()
            self.assertTrue(table.icm_rank.isna().all())
            self.assertEqual(int(table.f2000.max()), 1)
            coverage = reader.coverage()
            icm = coverage[coverage.group == 'icm']
            self.assertTrue(icm.numerator.isna().all())
            self.assertTrue(icm.rate.isna().all())
            self.assertEqual(driver.verify(config)['result'], 'PASS')

    def test_publish_preserves_the_retired_gene_axis_file_set(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            config.output.mkdir(parents=True, exist_ok=True)
            legacy = config.output / 'fused_rank.csv'
            legacy.write_text('retired product\n')
            driver.publish(config)
            self.assertTrue(legacy.is_file())

    def test_publish_is_deterministic_across_roots(self):
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            config_a = write_fixture(Path(first) / 'repo')
            config_b = write_fixture(Path(second) / 'repo')
            result_a = driver.publish(config_a)
            result_b = driver.publish(config_b)
            self.assertEqual(result_a['outputs'], result_b['outputs'])
            self.assertEqual(
                sha256_file(config_a.output / 'manifest.json'),
                sha256_file(config_b.output / 'manifest.json'))

    def test_verify_replays_and_detects_tampering(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            driver.publish(config)
            self.assertEqual(driver.verify(config)['result'], 'PASS')
            candidate = config.output / 'candidate_rank_sum.csv'
            original = candidate.read_text()
            candidate.write_text(original + '\n')
            with self.assertRaises(ValueError):
                driver.verify(config)
            candidate.write_text(original)

    def test_verify_refuses_a_pre_v2_release_with_a_self_describing_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            driver.publish(config)
            manifest_path = config.output / 'manifest.json'
            manifest = json.loads(manifest_path.read_text())
            manifest['release'].pop('algorithm')
            manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
            with self.assertRaises(ValueError) as raised:
                driver.verify(config)
            message = str(raised.exception)
            self.assertIn('predate the v2 f2000', message)
            self.assertIn('#208 must republish', message)

    def test_failed_publish_keeps_the_accepted_release(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            driver.publish(config)
            before = sha256_file(config.output / 'manifest.json')
            tenv = config.tenv_hvg
            original = tenv.read_text()
            corrupted = pd.read_csv(tenv)
            corrupted.loc[1, 'rank'] = corrupted.loc[0, 'rank']
            corrupted.to_csv(tenv, index=False, lineterminator='\n')
            with self.assertRaises(ValueError):
                driver.publish(config)
            self.assertEqual(sha256_file(config.output / 'manifest.json'), before)
            tenv.write_text(original)
            self.assertEqual(driver.verify(config)['result'], 'PASS')

    def test_uncontracted_fusion_shape_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            rewrite_fusion_output(config, 'nhf.deg2000.csv',
                                  lambda frame: frame.iloc[:-1])
            with self.assertRaisesRegex(ValueError, 'rows'):
                driver.publish(config)

def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()



class CoverageTests(unittest.TestCase):
    """Coverage and composition carry the DEG2000 numerator *and* denominator."""

    @staticmethod
    def inputs():
        axes = {4: np.array([0, 1, 6, 7]), 6: np.array([2, 3, 4, 5, 6, 7])}
        tenv = np.array([8, 7, 6, 5, 4, 3, 2, 1], dtype=np.int64)
        masks = [
            np.isin(np.arange(AXIS_SIZE), (0, 1)),   # nhf
            None,                                   # icm: excluded from the fusion
            np.isin(np.arange(AXIS_SIZE), (7,)),     # af: its single DEG2000 gene
            np.isin(np.arange(AXIS_SIZE), (2, 3, 4)),  # genetic_cardiomyopathy
            np.isin(np.arange(AXIS_SIZE), (0, 4)),   # cardiotox_TOP2
            np.isin(np.arange(AXIS_SIZE), (5,)),     # cardiotox_other
        ]
        sizes = [2, None, 1, 3, 2, 1]
        return axes, tenv, masks, sizes

    def test_coverage_carries_the_numerator_and_its_own_denominator(self):
        axes, tenv, masks, _sizes = self.inputs()
        table = source_topn_coverage(axes, tenv, masks, categories=CATEGORIES,
                                     sizes=(4, 6))
        self.assertEqual(list(table.columns), list(COVERAGE_COLUMNS))
        self.assertEqual(len(table), 2 * (1 + len(CATEGORIES)))
        rows = {record['group']: record for record in table.to_dict('records')
                if record['N'] == 4}
        # tenv-HVG keeps the top-N contract: |top_N(tenv) ∩ G_4| / N.
        self.assertEqual(int(rows['tenv_hvg']['numerator']), 2)
        self.assertEqual(int(rows['tenv_hvg']['denominator']), 4)
        self.assertEqual(float(rows['tenv_hvg']['rate']), 50.0)
        # A disease's denominator is its own DEG2000 size, not N and not the axis.
        self.assertEqual(int(rows['nhf']['numerator']), 2)
        self.assertEqual(int(rows['nhf']['denominator']), 2)
        self.assertEqual(float(rows['nhf']['rate']), 100.0)
        self.assertEqual(int(rows['genetic_cardiomyopathy']['numerator']), 0)
        self.assertEqual(int(rows['genetic_cardiomyopathy']['denominator']), 3)
        self.assertEqual(float(rows['genetic_cardiomyopathy']['rate']), 0.0)
        self.assertEqual(int(rows['cardiotox_TOP2']['numerator']), 1)
        self.assertEqual(int(rows['cardiotox_TOP2']['denominator']), 2)
        self.assertEqual(int(rows['cardiotox_other']['numerator']), 0)
        self.assertEqual(int(rows['cardiotox_other']['denominator']), 1)
        # AF's single DEG2000 gene 7 is on both fixture axes.
        self.assertEqual(int(rows['af']['numerator']), 1)
        self.assertEqual(int(rows['af']['denominator']), 1)
        self.assertEqual(float(rows['af']['rate']), 100.0)
        # An excluded disease keeps an empty row, never a zero.
        self.assertTrue(pd.isna(rows['icm']['numerator']))
        self.assertTrue(pd.isna(rows['icm']['denominator']))
        self.assertTrue(pd.isna(rows['icm']['rate']))
        larger = {record['group']: record for record in table.to_dict('records')
                  if record['N'] == 6}
        self.assertEqual(int(larger['nhf']['numerator']), 0)
        self.assertEqual(int(larger['genetic_cardiomyopathy']['numerator']), 3)
        self.assertEqual(float(larger['genetic_cardiomyopathy']['rate']), 100.0)
        for record in table.to_dict('records'):
            if pd.isna(record['numerator']):
                continue
            self.assertLessEqual(int(record['numerator']), int(record['denominator']))
            self.assertAlmostEqual(float(record['rate']),
                                   round(100.0 * int(record['numerator'])
                                         / int(record['denominator']), 2))

    def test_deg2000_composition_reports_the_numerator_and_denominator(self):
        _axes, _tenv, masks, sizes = self.inputs()
        composition = deg2000_composition(np.array([0, 1, 6, 7]), masks, sizes,
                                          categories=CATEGORIES)
        self.assertEqual(composition['nhf'],
                         {'n_covered': 2, 'n_deg2000': 2, 'share': 1.0})
        self.assertIsNone(composition['icm'])
        self.assertEqual(composition['af'],
                         {'n_covered': 1, 'n_deg2000': 1, 'share': 1.0})
        self.assertEqual(composition['genetic_cardiomyopathy'],
                         {'n_covered': 0, 'n_deg2000': 3, 'share': 0.0})
        self.assertEqual(composition['cardiotox_TOP2']['n_covered'], 1)
        self.assertEqual(composition['cardiotox_other']['n_covered'], 0)

    def test_f2000_distribution_is_bounded_by_the_active_disease_count(self):
        np.testing.assert_array_equal(f2000_distribution(np.array([0, 1, 1, 2]), 2),
                                      [1, 2, 1])
        np.testing.assert_array_equal(f2000_distribution(np.array([0, 0]),
                                                         n_diseases=1), [2, 0])
        with self.assertRaises(ValueError):
            f2000_distribution(np.array([0, 1, 2]), n_diseases=1)
        with self.assertRaises(ValueError):
            f2000_distribution(np.array([-1]))

    def test_out_of_contract_inputs_fail_closed(self):
        axes, tenv, masks, sizes = self.inputs()
        with self.assertRaises(ValueError):
            source_topn_coverage(axes, np.array([1, 1, 2, 3, 4, 5, 6, 7]),
                                 masks, categories=CATEGORIES, sizes=(4, 6))
        with self.assertRaises(ValueError):
            source_topn_coverage(axes, tenv, masks, categories=CATEGORIES, sizes=(5,))
        with self.assertRaises(ValueError):
            source_topn_coverage(axes, tenv, masks[:4], categories=CATEGORIES,
                                 sizes=(4,))
        with self.assertRaises(ValueError):
            source_topn_coverage(axes, tenv,
                                 [np.zeros(AXIS_SIZE, dtype=bool)] * 5,
                                 categories=CATEGORIES, sizes=(4,))
        with self.assertRaises(ValueError):
            source_topn_coverage(
                {4: np.array([0, 1, 6, 7])}, tenv, masks,
                categories=CATEGORIES, sizes=(4, 8))
        with self.assertRaises(ValueError):
            deg2000_composition(np.array([0, 1]), masks, [2, None, 1, 3],
                                categories=CATEGORIES)
        with self.assertRaises(ValueError):
            deg2000_composition(np.array([0, 1, 6, 7]),
                                [np.isin(np.arange(AXIS_SIZE), (0, 1)), None,
                                 np.isin(np.arange(AXIS_SIZE), (7,)),
                                 np.isin(np.arange(AXIS_SIZE), (2, 3, 4)),
                                 np.isin(np.arange(AXIS_SIZE), (0, 4)),
                                 np.isin(np.arange(AXIS_SIZE), (5,))],
                                [3, None, 1, 3, 2, 1], categories=CATEGORIES)

class FigureTests(unittest.TestCase):
    """The evidence PNG is a function of its plotted-values CSV alone."""

    @staticmethod
    def synthetic_rows():
        rows = [{'record': 'figure', 'title': 'fixture'},
                {'record': 'panel', 'panel': 'a', 'axis': '4',
                 'panel_title': 'N = 4'},
                {'record': 'panel', 'panel': 'b', 'axis': '6',
                 'panel_title': 'N = 6'}]
        counts = (10, 11, 12, 13)
        for panel, selected in (('a', (1, 1, 1, 1)), ('b', (2, 2, 1, 1))):
            for index, (count, chosen) in enumerate(zip(counts, selected)):
                rows.append({'record': 'hist', 'panel': panel,
                             'bin_left': index, 'bin_right': index + 1,
                             'bin_center': index + 0.5, 'count': count,
                             'selected': chosen})
        return rows

    def test_renderer_rejects_incomplete_csv(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'bad.csv'
            path.write_text('record,panel,panel_title,title\nfigure,,,,x\n')
            with self.assertRaises(ValueError):
                renderer.load(path)

    def test_renderer_rejects_impossible_selection(self):
        def write(rows, path):
            fields = sorted({key for row in rows for key in row})
            with path.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)

        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            over = self.synthetic_rows()
            over[-1]['selected'] = over[-1]['count'] + 1
            write(over, temporary / 'over.csv')
            with self.assertRaises(ValueError):
                renderer.load(temporary / 'over.csv')
            short = self.synthetic_rows()
            for row in short:
                if row.get('panel') == 'a' and row.get('record') == 'hist':
                    row['selected'] = 0
            write(short, temporary / 'short.csv')
            with self.assertRaises(ValueError):
                renderer.load(temporary / 'short.csv')

    def test_synthetic_csv_renders_byte_identically_twice(self):
        rows = self.synthetic_rows()
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            path = temporary / 'figure.csv'
            fields = sorted({key for row in rows for key in row})
            with path.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            first, second = temporary / 'first.png', temporary / 'second.png'
            renderer.render(*renderer.load(path), first)
            renderer.render(*renderer.load(path), second)
            self.assertEqual(sha256_file(first), sha256_file(second))

    def test_committed_figure_replays_from_its_csv_alone(self):
        csv_path = STEP / 'figure/data/gene_axis_rank_sum_histogram.csv'
        png_path = STEP / 'figure/gene_axis_rank_sum_histogram.png'
        self.assertTrue(csv_path.is_file(), f'missing plotted-values CSV: {csv_path}')
        self.assertTrue(png_path.is_file(), f'missing committed figure: {png_path}')
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            copied = temporary / 'figure.csv'
            copied.write_bytes(csv_path.read_bytes())
            rendered = temporary / 'rendered.png'
            renderer.render(*renderer.load(copied), rendered)
            self.assertEqual(sha256_file(rendered), sha256_file(png_path))
        source = Path(renderer.__file__).read_text(encoding='utf-8')
        self.assertNotIn('read_csv', source)
        self.assertNotIn('read_excel', source)
        self.assertNotIn("palette.csv", source)
        self.assertNotIn("manifest.json", source)


class SourceCoverageFigureTests(unittest.TestCase):
    """The two-panel source-coverage PNG is a function of its plotted-values CSV alone."""

    @staticmethod
    def synthetic_rows():
        rows = [{'record': 'panel', 'panel': 'a',
                 'panel_title': 'Covered Gene Count', 'y_label': 'covered genes',
                 'value': 'count'},
                {'record': 'panel', 'panel': 'b',
                 'panel_title': 'Top-N Coverage Rate',
                 'y_label': 'coverage rate (%)', 'value': 'rate'}]
        for series, label, counts, rates in (
                ('G_4', 'G_4 (N = 4)', (2, 1), (50.0, 25.0)),
                ('G_6', 'G_6 (N = 6)', (6, 4), (100.0, 66.67))):
            for source, source_label, count, rate in zip(
                    ('tenv_hvg', 'nhf'), ('tenv-HVG', 'nHF'), counts, rates):
                rows.append({'record': 'bar', 'series': series,
                             'series_label': label, 'source': source,
                             'source_label': source_label,
                             'count': str(count), 'rate': str(rate)})
        return rows

    @staticmethod
    def write(rows, path):
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with Path(path).open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    def test_renderer_rejects_incomplete_or_impossible_csv(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            empty = temporary / 'empty.csv'
            empty.write_text('record,panel,panel_title,y_label,value,series,'
                             'series_label,source,source_label,count,rate\n')
            with self.assertRaises(ValueError):
                coverage_renderer.load(empty)
            # Bars without the panel records that say what to draw.
            self.write(self.synthetic_rows()[2:], temporary / 'panels.csv')
            with self.assertRaises(ValueError):
                coverage_renderer.load(temporary / 'panels.csv')
            self.write(self.synthetic_rows()[:-1], temporary / 'missing.csv')
            with self.assertRaises(ValueError):
                coverage_renderer.load(temporary / 'missing.csv')
            # A panel that names a column the bars do not carry.
            rows = self.synthetic_rows()
            rows[0]['value'] = 'unknown'
            self.write(rows, temporary / 'column.csv')
            with self.assertRaises(ValueError):
                coverage_renderer.load(temporary / 'column.csv')
            for name, key, value in (('negative', 'count', '-1'),
                                     ('fractional', 'count', '1.5'),
                                     ('rate', 'rate', '101')):
                rows = self.synthetic_rows()
                rows[2][key] = value
                self.write(rows, temporary / f'{name}.csv')
                with self.assertRaises(ValueError):
                    coverage_renderer.load(temporary / f'{name}.csv')

    def test_synthetic_csv_renders_byte_identically_twice(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            path = temporary / 'figure.csv'
            self.write(self.synthetic_rows(), path)
            first, second = temporary / 'first.png', temporary / 'second.png'
            coverage_renderer.render(*coverage_renderer.load(path), first)
            coverage_renderer.render(*coverage_renderer.load(path), second)
            self.assertEqual(sha256_file(first), sha256_file(second))

    def test_committed_figure_replays_from_its_csv_alone(self):
        csv_path = STEP / 'figure/data/gene_axis_source_topn_coverage.csv'
        png_path = STEP / 'figure/gene_axis_source_topn_coverage.png'
        self.assertTrue(csv_path.is_file(), f'missing plotted-values CSV: {csv_path}')
        self.assertTrue(png_path.is_file(), f'missing committed figure: {png_path}')
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / 'figure.csv'
            copied.write_bytes(csv_path.read_bytes())
            rendered = Path(temporary) / 'rendered.png'
            coverage_renderer.render(*coverage_renderer.load(copied), rendered)
            self.assertEqual(sha256_file(rendered), sha256_file(png_path))
        source = Path(coverage_renderer.__file__).read_text(encoding='utf-8')
        for token in ('read_csv', 'read_excel', 'palette.csv', 'manifest.json'):
            self.assertNotIn(token, source)




def source_table_row(document, catalogs, dataset_id):
    """Return the cells of ``dataset_id``'s source-table row (row order = catalog order)."""
    offset = 0
    for category in gene_axis.RELEASE_CATEGORIES:
        catalog = json.loads(Path(catalogs[category]).read_text(encoding='utf-8'))
        for dataset in catalog['datasets']:
            if dataset['dataset_id'] == dataset_id:
                return [cell.text for cell in document.tables[0].rows[1 + offset].cells]
            offset += 1
    raise AssertionError(f'{dataset_id} is not in the concatenated catalogs')


class DocxTests(unittest.TestCase):
    """The rule document replays byte-identically and every table cell is contracted."""

    def catalogs(self):
        return {category: STEP / 'output/disease-deg' / category / 'catalog.json'
                for category in gene_axis.RELEASE_CATEGORIES}

    def test_committed_docx_replays_from_its_inputs(self):
        coverage = STEP / 'output/gene-axis/source_topn_coverage.csv'
        figure = STEP / 'figure/gene_axis_rank_sum_histogram.png'
        fusion = STEP / 'output/disease-fusion/catalog.json'
        axis = STEP / 'output/gene-identity/genes.csv'
        manifest = STEP / 'output/gene-axis/manifest.json'
        catalogs = self.catalogs()
        committed = STEP / 'output/gene-axis/gene_axis_fill_rule.docx'
        for path in (coverage, figure, fusion, axis, manifest, committed,
                     *catalogs.values()):
            self.assertTrue(path.is_file(), f'missing docx input or output: {path}')
        import gene_axis_source_metadata as metadata
        import make_gene_axis_docx as docx_builder
        with tempfile.TemporaryDirectory() as temporary:
            rendered = Path(temporary) / 'gene_axis_fill_rule.docx'
            docx_builder.build(coverage, figure, rendered, catalogs=catalogs,
                               metadata=metadata, fusion_catalog=fusion,
                               axis_path=axis, manifest_path=manifest)
            self.assertEqual(sha256_file(rendered), sha256_file(committed))

    def test_build_reads_the_five_catalog_paths_it_is_given(self):
        """A changed catalog ``n_deg`` must reach the rendered cell."""
        from docx import Document
        import make_gene_axis_docx as docx_builder
        coverage = STEP / 'output/gene-axis/source_topn_coverage.csv'
        figure = STEP / 'figure/gene_axis_rank_sum_histogram.png'
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            catalogs = {}
            for category, path in self.catalogs().items():
                copied = temporary / category / 'catalog.json'
                copied.parent.mkdir(parents=True)
                shutil.copyfile(path, copied)
                catalogs[category] = copied
            catalog = json.loads(catalogs['cardiotox_TOP2'].read_text(encoding='utf-8'))
            for dataset in catalog['datasets']:
                if dataset['dataset_id'] == 'gse243674':
                    dataset['n_deg'] = 12345
            catalogs['cardiotox_TOP2'].write_text(
                json.dumps(catalog, indent=2, sort_keys=True) + '\n', encoding='utf-8')
            rendered = temporary / 'wired.docx'
            docx_builder.build(coverage, figure, rendered, catalogs=catalogs)
            document = Document(str(rendered))
            cells = source_table_row(document, catalogs, 'gse243674')
            self.assertEqual(cells[6].replace(',', ''), '12345',
                             'the build ignored the catalog path it was handed')

    def test_build_reads_the_metadata_module_it_is_given(self):
        """A changed technology string must reach the rendered cell."""
        from docx import Document
        import make_gene_axis_docx as docx_builder
        module_path = EXECUTION / 'gene_axis_source_metadata.py'
        coverage = STEP / 'output/gene-axis/source_topn_coverage.csv'
        figure = STEP / 'figure/gene_axis_rank_sum_histogram.png'
        tracked = "'bulk RNA-seq (hiPSC-CM monolayer, 3 h / 24 h)'"
        probe = "'bulk RNA-seq (hiPSC-CM monolayer, 3 h / 24 h; metadata probe)'"
        source = module_path.read_text(encoding='utf-8')
        self.assertEqual(source.count(tracked), 1, 'fixture assumption: one matching technology')
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            copied = temporary / 'gene_axis_source_metadata.py'
            copied.write_text(source.replace(tracked, probe), encoding='utf-8')
            metadata = docx_builder.load_metadata_module(copied)
            rendered = temporary / 'wired.docx'
            docx_builder.build(coverage, figure, rendered,
                               catalogs=self.catalogs(), metadata=metadata)
            document = Document(str(rendered))
            cells = source_table_row(document, self.catalogs(), 'gse243674')
            self.assertEqual(cells[5], probe.strip("'"),
                             'the build ignored the metadata module it was handed')

    def test_cli_catalog_and_metadata_overrides_are_read(self):
        """The ``--catalog``/``--metadata`` CLI inputs reach the document too."""
        from docx import Document
        coverage = STEP / 'output/gene-axis/source_topn_coverage.csv'
        figure = STEP / 'figure/gene_axis_rank_sum_histogram.png'
        tracked = "'bulk RNA-seq (hiPSC-CM monolayer, 3 h / 24 h)'"
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            catalog_path = temporary / 'nhf-catalog.json'
            shutil.copyfile(self.catalogs()['nhf'], catalog_path)
            catalog = json.loads(catalog_path.read_text(encoding='utf-8'))
            catalog['datasets'][0]['n_deg'] = 424242
            catalog_path.write_text(
                json.dumps(catalog, indent=2, sort_keys=True) + '\n', encoding='utf-8')
            metadata_path = temporary / 'gene_axis_source_metadata.py'
            metadata_path.write_text(
                (EXECUTION / 'gene_axis_source_metadata.py')
                .read_text(encoding='utf-8')
                .replace(tracked, tracked[:-2] + "; cli probe)'"),
                encoding='utf-8')
            out = temporary / 'cli.docx'
            result = subprocess.run(
                [sys.executable, str(EXECUTION / 'make_gene_axis_docx.py'),
                 '--coverage', str(coverage), '--figure', str(figure),
                 '--catalog', f'nhf={catalog_path}',
                 '--metadata', str(metadata_path), '--out', str(out)],
                capture_output=True, text=True, timeout=300)
            self.assertEqual(result.returncode, 0, result.stderr)
            document = Document(str(out))
            first_dataset = json.loads(
                self.catalogs()['nhf'].read_text(encoding='utf-8'))['datasets'][0]
            cells = source_table_row(document, self.catalogs(), first_dataset['dataset_id'])
            self.assertEqual(cells[6].replace(',', ''), '424242', '--catalog was ignored')
            cells = source_table_row(document, self.catalogs(), 'gse243674')
            self.assertEqual(cells[5], tracked.strip("'")[:-1] + '; cli probe)',
                             '--metadata was ignored')

    def test_coverage_table_columns_are_contracted(self):
        """The long-format columns plus the v2 denominators, group by group."""
        coverage = STEP / 'output/gene-axis/source_topn_coverage.csv'
        import make_gene_axis_docx as docx_builder
        table = docx_builder.read_coverage(coverage)
        self.assertEqual(list(table.columns), list(docx_builder.COVERAGE_COLUMNS))
        fusion = json.loads((STEP / 'output/disease-fusion/catalog.json').read_text(
            encoding='utf-8'))
        deg2000 = {record['category']: record['deg2000'] for record in fusion['diseases']}
        for record in table.itertuples(index=False):
            with self.subTest(N=record.N, group=record.group):
                if record.group == 'tenv_hvg':
                    self.assertEqual(record.denominator, record.N,
                                     'tenv-HVG divides by its top-N set size')
                else:
                    label = gene_axis.CATEGORY_LABELS[record.group]
                    self.assertEqual(record.denominator, deg2000[label],
                                     'a disease divides by its own |DEG2000_d|')
                self.assertEqual(record.n_axis, record.N, 'n_axis is the released axis size')

    def test_source_table_rows_match_catalogs(self):
        """The first table is the released source dataset table, row for row.

        The v2 release publishes 27 dataset rows across six categories; the
        tenv-HVG fusion row makes the 28th.
        """
        from docx import Document
        import gene_axis_source_metadata as metadata
        import make_gene_axis_docx as docx_builder
        committed = STEP / 'output/gene-axis/gene_axis_fill_rule.docx'
        self.assertTrue(committed.is_file(), f'missing committed docx: {committed}')
        document = Document(str(committed))
        self.assertGreaterEqual(len(document.tables), 1, 'the document has no tables')
        table = document.tables[0]
        header = [cell.text for cell in table.rows[0].cells]
        self.assertEqual(header, list(docx_builder.SOURCE_TABLE_COLUMNS))
        self.assertEqual(len(table.rows), 1 + 28, 'expected 28 source rows + header')
        expected = []
        for category in gene_axis.RELEASE_CATEGORIES:
            catalog = json.loads(self.catalogs()[category].read_text(encoding='utf-8'))
            expected.extend((catalog['category'], dataset) for dataset in catalog['datasets'])
        self.assertEqual(len(expected), 27)
        for offset, (category_label, dataset) in enumerate(expected, start=1):
            cells = [cell.text for cell in table.rows[offset].cells]
            dataset_id = dataset['dataset_id']
            source = metadata.SOURCES[dataset_id]
            with self.subTest(dataset=dataset_id):
                self.assertEqual(cells[0], category_label, 'category column')
                self.assertEqual(cells[2], dataset['model_system'], 'sample type column')
                groups = dataset.get('arm_or_class')
                if groups is not None:
                    self.assertEqual(cells[3], groups, 'groups column')
                else:
                    # Cardiotox rows: class label/id/drugs plus the vehicle phrase the
                    # catalog's own ``shared_control_group`` publishes.
                    shared = dataset['shared_control_group']
                    _, separator, tail = shared.partition(':')
                    self.assertTrue(separator, f'{dataset_id}: no ":" in shared_control_group')
                    vehicle = tail.strip().split(' (', 1)[0]
                    expected_groups = (
                        f"{dataset['class_label']} ({dataset['class_id']}; "
                        f"{', '.join(dataset['member_drugs'])}) vs {vehicle}")
                    self.assertEqual(cells[3], expected_groups,
                                     'groups cell must derive from the catalog')
                self.assertEqual(
                    cells[4], f"{dataset['n_case']} / {dataset['n_control']}",
                    'sample sizes column')
                self.assertEqual(cells[5], source['technology'], 'sequencing technology column')
                accession = source['accession']
                expected_accession = (f"{accession} ({source['accession_note']})"
                                      if source['accession_note'] else accession)
                self.assertEqual(cells[1], expected_accession, 'accession column')
                self.assertEqual(int(cells[6].replace(',', '')), dataset['n_deg'],
                                 'CM DEG count column is the v2 DEG count, not n_significant')
        last = [cell.text for cell in table.rows[-1].cells]
        tenv = metadata.SOURCES[metadata.TENV_HVG]
        self.assertEqual(last[0], 'tenv-HVG (fusion source)')
        self.assertLessEqual(len(last[1]), 60, 'tenv accession cell must stay short')
        self.assertIn('GSE262419', last[1])
        self.assertIn(tenv['technology'], last[5])
        self.assertIn('not a DEG source', last[6])

    def test_sample_size_footnote_cross_checks_the_mapping(self):
        """The footnote's five sample-size bases must exist in the tracked notes."""
        from docx import Document
        import gene_axis_source_metadata as metadata
        committed = STEP / 'output/gene-axis/gene_axis_fill_rule.docx'
        self.assertTrue(committed.is_file(), f'missing committed docx: {committed}')
        document = Document(str(committed))
        footnote = next((paragraph.text for paragraph in document.paragraphs
                         if '样本量口径' in paragraph.text), None)
        self.assertIsNotNone(footnote, 'the source-table footnote lacks the sample-size basis')
        groups = {
            '46': ('gse270788__nicm', 'gse270788__icm', 'gse270788__ami'),
            '1,171': ('dtoxs_ipsc', 'dtoxs_ipsc__atc_l4-L01EA',
                      'dtoxs_ipsc__atc_l4-L01EB', 'dtoxs_ipsc__atc_l4-L01EX'),
            '72': ('gse243674',),
            '3,072': ('gse244740',),
            '79': ('reichart-cm__dcm_genotype_negative',
                   'reichart-cm__dcm_genotype_positive', 'reichart-cm__arvc'),
            '17': ('gse35108__dcm_ipsc_cm',),
            '177': ('gse165303__dcm_lv_bulk',),
            '48': ('gse138252__af',),
            '14': ('gse238242__af',),
            '191,795': ('kuppe-mi__iz',),
        }
        for number, dataset_ids in groups.items():
            with self.subTest(number=number):
                self.assertIn(number, footnote, f'footnote lacks the {number} basis')
                notes = [metadata.SOURCES[dataset_id]['sample_size_note']
                         for dataset_id in dataset_ids]
                self.assertTrue(
                    any(note and number in note for note in notes),
                    f'the mapping notes no longer record {number}')
                for dataset_id in dataset_ids:
                    note = metadata.SOURCES[dataset_id]['sample_size_note']
                    if note is not None:
                        self.assertIsInstance(note, str)
                        self.assertTrue(note.strip(),
                                        f'{dataset_id}: blank sample_size_note')

    def test_disease_table_matches_the_published_catalogs(self):
        """The second table is the per-disease DEG/DEG2000 record, row for row."""
        from docx import Document
        import make_gene_axis_docx as docx_builder
        committed = STEP / 'output/gene-axis/gene_axis_fill_rule.docx'
        self.assertTrue(committed.is_file(), f'missing committed docx: {committed}')
        document = Document(str(committed))
        self.assertGreaterEqual(len(document.tables), 2,
                                'the document has no per-disease table')
        table = document.tables[1]
        header = [cell.text for cell in table.rows[0].cells]
        self.assertEqual(header, list(docx_builder.DISEASE_TABLE_COLUMNS))
        self.assertEqual(len(table.rows), 1 + 6, 'expected 6 disease rows + header')
        fusion = json.loads((STEP / 'output/disease-fusion/catalog.json').read_text(
            encoding='utf-8'))
        by_label = {record['category']: record for record in fusion['diseases']}
        rows = {row.cells[0].text: [cell.text for cell in row.cells]
                for row in table.rows[1:]}
        self.assertEqual(list(rows), [gene_axis.CATEGORY_LABELS[category]
                                      for category in gene_axis.RELEASE_CATEGORIES])
        for category in gene_axis.RELEASE_CATEGORIES:
            label = gene_axis.CATEGORY_LABELS[category]
            catalog = json.loads(self.catalogs()[category].read_text(encoding='utf-8'))
            record = by_label[label]
            cells = rows[label]
            with self.subTest(disease=label):
                self.assertEqual(int(cells[1]), record['declared_source_count'],
                                 'source-dataset count')
                self.assertEqual(int(cells[3].replace(',', '')), record['n_deg'],
                                 'DEG union')
                self.assertEqual(int(cells[4].replace(',', '')), record['deg2000'],
                                 'DEG2000')
                self.assertEqual(int(cells[4].replace(',', '')), catalog['deg2000'],
                                 'the two published catalogs must agree')
                tiers = sorted({entry['tier'] for entry in catalog['tier_sources']})
                self.assertEqual(cells[2], ' → '.join(str(tier) for tier in tiers),
                                 'tier path')
        # The uncomfortable results are stated, never smoothed over.
        self.assertIn('SLC9A7', rows['AF'][5], 'AF must name its single gene')
        self.assertIn('未达到', rows['AF'][5])
        self.assertIn('唯一达到', rows['cardiotox_TOP2'][5])
        for short in ('nHF', 'ICM'):
            self.assertIn('未达到', rows[short][5], f'{short} must state its shortfall')
            self.assertIn('没有更低层可收', rows[short][5])

    def test_rule_section_opens_the_document_with_the_v2_rule(self):
        """A reader meets the v2 rule before any table, and never the pN rule."""
        from docx import Document
        import make_gene_axis_docx as docx_builder
        committed = STEP / 'output/gene-axis/gene_axis_fill_rule.docx'
        document = Document(str(committed))
        headings = [paragraph.text for paragraph in document.paragraphs
                    if paragraph.style.name == 'Title'
                    or paragraph.style.name.startswith('Heading')]
        self.assertEqual(headings[0], docx_builder.TITLE)
        self.assertEqual(headings[1], docx_builder.RULE_HEADING,
                         'the v2 rule is the opening section')
        text = '\n'.join(paragraph.text for paragraph in document.paragraphs)
        for token in ('padj < 0.05', '|log2FC| ≥ 2', 'pass-any-contrast', 'DEG2000',
                      'f2000', 'tenv-HVG', 'tier 1', 'tier 3'):
            self.assertIn(token, text)
        self.assertIn('pN 已退役', text)
        self.assertNotIn('按 pN 降序', text, 'the retired selection must not be described')


if __name__ == '__main__':
    unittest.main()
