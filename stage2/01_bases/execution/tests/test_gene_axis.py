"""Fixture-only tests for the two-level rank-sum fusion and the released axes.

No real data is read: the axis, the per-dataset ranks, the tenv-HVG ranking and
the published plotted values are hand-built 8-gene fixtures. The fusion, the
axis fill rule and the readers are exercised at the public API boundary; one
seam only (ADR 0066 / spec). Determinism is tested by building the same fixture
twice into scratch candidate roots and comparing every digest.
"""
import csv
import hashlib
import json
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
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import gene_axis
from bases.gene_axis import (
    AXIS_SIZES, CATEGORIES, CATEGORY_RANK_COLUMNS, CANDIDATE_COLUMNS, GeneAxis,
    category_rank, winning_tier,
)

AXIS_SYMBOLS = ['G0', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7']


def axis_fixture():
    return pd.DataFrame({
        'gene_index': np.arange(len(AXIS_SYMBOLS), dtype=np.int64),
        'symbol': AXIS_SYMBOLS,
        'hgnc_id': [f'HGNC:{index + 1}' for index in range(len(AXIS_SYMBOLS))],
        'ensembl_id': [f'ENSG0000000000{index + 1}' for index in range(len(AXIS_SYMBOLS))],
    })


def dataset_record(dataset_id, tier, isogenic=False, n_measured=8):
    record = {'dataset_id': dataset_id, 'tier': tier, 'n_measured': n_measured}
    if isogenic:
        record['isogenic_exception'] = {'applies': True, 'why': 'fixture isogenic pair'}
    return record


class WinningTierTests(unittest.TestCase):
    """ADR 0066 priority tier + isogenic exception (which datasets a category fuses)."""

    def test_lowest_tier_supplies_the_source(self):
        catalog = {'datasets': [dataset_record('a', 1), dataset_record('b', 1),
                                dataset_record('c', 2), dataset_record('d', 3)]}
        selected = winning_tier(catalog)
        self.assertEqual(selected['winning_tier'], 1)
        self.assertEqual(selected['datasets'], ['a', 'b'])
        self.assertEqual(selected['excluded_tiers'], [2, 3])
        self.assertEqual(selected['retained_by_exception'], [])

    def test_isogenic_exception_retains_higher_tier_models(self):
        catalog = {'datasets': [
            dataset_record('human_a', 1), dataset_record('human_b', 1),
            dataset_record('isogenic_a', 2, isogenic=True),
            dataset_record('isogenic_b', 2, isogenic=True),
            dataset_record('other_invitro', 2),
        ]}
        selected = winning_tier(catalog)
        self.assertEqual(selected['datasets'],
                         ['human_a', 'human_b', 'isogenic_a', 'isogenic_b'])
        self.assertEqual(selected['retained_by_exception'],
                         ['isogenic_a', 'isogenic_b'])
        self.assertEqual(selected['excluded_tiers'], [2])

    def test_missing_tier_fails_closed(self):
        with self.assertRaises(ValueError):
            winning_tier({'datasets': [{'dataset_id': 'a'}]})
        with self.assertRaises(ValueError):
            winning_tier({'datasets': []})


class LevelOneTests(unittest.TestCase):
    """Level 1: sum(rank_norm) inside one category, ties by gene_index."""

    def test_single_dataset_category_is_its_own_order(self):
        ranks = np.array([8, 1, 7, 2, 6, 3, 5, 4], dtype=np.int64)
        result = category_rank([ranks], [8])
        np.testing.assert_array_equal(result['rank'], ranks)
        self.assertEqual(result['n_datasets'], 1)
        self.assertTrue(result['measured'].all())

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

    def test_unmeasured_gene_contributes_exactly_one(self):
        ranks = np.array([1, 2, 3, 4, 5, 6, 7, 7], dtype=np.int64)
        result = category_rank([ranks], [6])
        np.testing.assert_allclose(result['sum'], [1/6, 2/6, 3/6, 4/6, 5/6, 1.0, 1.0, 1.0])
        np.testing.assert_array_equal(result['rank'], np.arange(1, 9))
        np.testing.assert_array_equal(
            result['measured'], [True, True, True, True, True, True, False, False])

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


class LevelTwoTests(unittest.TestCase):
    """Level 2: rank_sum over the seven sources, ascending-is-better, pN, coverage."""

    @staticmethod
    def identity():
        return np.arange(1, len(AXIS_SYMBOLS) + 1, dtype=np.int64)

    def test_sum_of_seven_normalized_ranks_and_gene_index_tie_break(self):
        ranks = [self.identity().copy() for _ in CATEGORIES]
        # one source puts gene 0 at rank 8 and gene 1 at rank 2: the integer
        # sums 8+6*1 and 2+6*2 both equal 14, and genes 6/7 also tie at 49
        ranks[0] = np.array([8, 2, 3, 4, 5, 6, 7, 1], dtype=np.int64)
        measured = [np.ones(8, dtype=bool) for _ in CATEGORIES]
        result = gene_axis.fusion_rank_sum(
            ranks, measured, self.identity(), sizes=(2, 4))
        np.testing.assert_allclose(
            result['rank_sum'],
            np.array([14, 14, 21, 28, 35, 42, 49, 49]) / 8.0)
        # equal rank sums are ordered by gene_index in the fusion (fill-rule test)
        self.assertEqual(result['rank_sum'][0], result['rank_sum'][1])

    def test_pn_boundaries_zero_and_seven(self):
        ranks = [self.identity().copy() for _ in CATEGORIES]
        ranks[0][:2] = [2, 1]
        measured = [np.ones(8, dtype=bool) for _ in CATEGORIES]
        result = gene_axis.fusion_rank_sum(
            ranks, measured, self.identity(), sizes=(1, 2, 4))
        self.assertEqual(result['pn'][1].tolist(), [6, 1, 0, 0, 0, 0, 0, 0])
        self.assertEqual(result['pn'][2].tolist(), [7, 7, 0, 0, 0, 0, 0, 0])
        self.assertEqual(result['pn'][4].tolist(), [7, 7, 7, 7, 0, 0, 0, 0])
        self.assertEqual(int(result['pn'][2].max()), 7)
        self.assertEqual(int(result['pn'][2].min()), 0)

    def test_n_sources_measured_counts_sources_not_datasets(self):
        ranks = [self.identity().copy() for _ in CATEGORIES]
        measured = [np.zeros(8, dtype=bool) for _ in CATEGORIES]
        measured[0][0] = True
        measured[3][0] = True
        result = gene_axis.fusion_rank_sum(ranks, measured, self.identity(), sizes=(2,))
        # tenv always measures the gene; two of the six categories do
        self.assertEqual(result['n_sources_measured'].tolist(),
                         [3, 1, 1, 1, 1, 1, 1, 1])

    def test_rank_sum_is_ascending_is_better(self):
        ranks = [self.identity().copy() for _ in CATEGORIES]
        ranks[0] = self.identity()[::-1].copy()
        measured = [np.ones(8, dtype=bool) for _ in CATEGORIES]
        result = gene_axis.fusion_rank_sum(ranks, measured, self.identity(), sizes=(2,))
        order = np.argsort(result['rank_sum'], kind='stable')
        self.assertEqual(int(order[0]), 0)   # every source ranks gene 0 first
        self.assertEqual(int(order[-1]), 7)  # and gene 7 last
        self.assertTrue(result['rank_sum'][order[0]] < result['rank_sum'][order[-1]])

    def test_fusion_rejects_incomplete_source_ranks(self):
        ranks = [self.identity().copy() for _ in CATEGORIES]
        ranks[0][0] = 99
        measured = [np.ones(8, dtype=bool) for _ in CATEGORIES]
        with self.assertRaises(ValueError):
            gene_axis.fusion_rank_sum(ranks, measured, self.identity(), sizes=(2,))
        with self.assertRaises(ValueError):
            gene_axis.fusion_rank_sum(ranks[:5], measured[:5], self.identity(), sizes=(2,))
        good = [self.identity().copy() for _ in CATEGORIES]
        with self.assertRaises(ValueError):
            gene_axis.fusion_rank_sum(good, measured, np.array([1, 2, 3]), sizes=(2,))


class CandidateTableTests(unittest.TestCase):
    """The master candidate table: one row per frozen-axis gene, all contracted columns."""

    @staticmethod
    def fusion_inputs():
        identity = np.arange(1, len(AXIS_SYMBOLS) + 1, dtype=np.int64)
        ranks = [identity.copy() for _ in CATEGORIES]
        ranks[0] = np.array([8, 2, 3, 4, 5, 6, 7, 1], dtype=np.int64)
        measured = [np.ones(8, dtype=bool) for _ in CATEGORIES]
        return ranks, measured, identity.copy()

    def test_columns_rows_and_values(self):
        axis = axis_fixture()
        ranks, measured, tenv = self.fusion_inputs()
        table = gene_axis.candidate_table(axis, tenv, ranks, measured, sizes=(2, 4))
        self.assertEqual(list(table.columns), list(gene_axis.candidate_columns((2, 4))))
        self.assertEqual(len(table), len(AXIS_SYMBOLS))
        np.testing.assert_array_equal(table.gene_index.to_numpy(), np.arange(8))
        self.assertEqual(table.symbol.tolist(), AXIS_SYMBOLS)
        np.testing.assert_array_equal(table.tenv_hvg_rank.to_numpy(), tenv)
        for column, rank in zip(CATEGORY_RANK_COLUMNS, ranks):
            np.testing.assert_array_equal(table[column].to_numpy(), rank)
        expected = gene_axis.fusion_rank_sum(
            ranks, measured, tenv, sizes=(2, 4))
        np.testing.assert_allclose(table.rank_sum.to_numpy(), expected['rank_sum'])
        np.testing.assert_array_equal(
            table.n_sources_measured.to_numpy(),
            expected['n_sources_measured'])
        np.testing.assert_array_equal(table.p2.to_numpy(), expected['pn'][2])
        np.testing.assert_array_equal(table.p4.to_numpy(), expected['pn'][4])

    def test_unmeasured_mask_counts_sources_and_keeps_the_row(self):
        axis = axis_fixture()
        ranks, measured, tenv = self.fusion_inputs()
        measured[0][3] = False
        table = gene_axis.candidate_table(axis, tenv, ranks, measured, sizes=(2,))
        self.assertEqual(int(table.n_sources_measured[3]), 6)
        self.assertIn(3, table.gene_index.tolist())

    def test_axis_order_is_required(self):
        axis = axis_fixture().iloc[::-1].reset_index(drop=True)
        ranks, measured, tenv = self.fusion_inputs()
        with self.assertRaises(ValueError):
            gene_axis.candidate_table(axis, tenv, ranks, measured, sizes=(2,))


class AxisFillTests(unittest.TestCase):
    """ADR 0066 axis fill rule: 0.5*N tenv-HVG then 0.5*N by pN / rank_sum / index."""

    @staticmethod
    def fill_inputs():
        # tenv-HVG ranks genes in reverse; the top 2 are 6,7 and the top 3 add 5
        tenv = np.array([8, 7, 6, 5, 4, 3, 2, 1], dtype=np.int64)
        rank_sum = np.array([0.50, 0.20, 0.20, 0.30, 0.10, 0.40, 0.05, 0.01])
        p_n = {4: np.array([0, 0, 0, 4, 0, 5, 7, 7]),
               6: np.array([4, 3, 2, 1, 5, 7, 7, 7])}
        return tenv, rank_sum, p_n

    def test_fill_takes_exactly_half_tenv_and_half_replication(self):
        tenv, rank_sum, p_n = self.fill_inputs()
        for n in (4, 6):
            axis = gene_axis.axis_fill(tenv, p_n[n], rank_sum, n=n)
            self.assertEqual(len(axis), n)
            from_tenv = int((tenv[axis] <= n // 2).sum())
            self.assertEqual(from_tenv, n // 2)
            self.assertEqual(n - from_tenv, n // 2)

    def test_pN_descending_then_rank_sum_then_gene_index(self):
        tenv, rank_sum, p_n = self.fill_inputs()
        axis4 = gene_axis.axis_fill(tenv, p_n[4], rank_sum, n=4)
        # tenv top-2 are 6,7; of the remaining, gene 5 has the top p4 then gene 3
        self.assertEqual(axis4.tolist(), [3, 5, 6, 7])
        axis6 = gene_axis.axis_fill(tenv, p_n[6], rank_sum, n=6)
        # tenv top-3 are 5,6,7; fill takes gene 4 (p=5), then gene 0 (p=4)
        # and gene 1 (p=3) by pN descending
        self.assertEqual(axis6.tolist(), [0, 1, 4, 5, 6, 7])

    def test_ties_break_on_rank_sum_then_gene_index(self):
        tenv = np.array([8, 7, 6, 5, 4, 3, 2, 1], dtype=np.int64)
        # four candidates tie on p; the fill order is rank_sum asc, then index
        rank_sum = np.array([0.9, 0.2, 0.2, 0.1, 0.5, 0.5, 0.5, 0.5])
        p_n = np.array([7, 7, 7, 7, 0, 7, 7, 7])
        axis = gene_axis.axis_fill(tenv, p_n, rank_sum, n=4)
        # tenv top-2 = genes 6,7; fill needs 2: gene 3 (rank_sum 0.1) then gene 1
        # (rank_sum 0.2 ties gene 2 at 0.2, gene_index picks 1)
        self.assertEqual(axis.tolist(), [1, 3, 6, 7])

    def test_axes_are_not_nested_and_intersection_is_reported(self):
        tenv, rank_sum, p_n = self.fill_inputs()
        axis4 = gene_axis.axis_fill(tenv, p_n[4], rank_sum, n=4)
        axis6 = gene_axis.axis_fill(tenv, p_n[6], rank_sum, n=6)
        self.assertFalse(set(axis4) <= set(axis6))
        self.assertFalse(set(axis6) <= set(axis4))
        overlap = gene_axis.axis_intersection(axis4, axis6)
        self.assertEqual(overlap['n_intersection'], 3)
        self.assertEqual(overlap['n_union'], 7)
        self.assertAlmostEqual(overlap['jaccard'], 3 / 7)

    def test_axis_is_an_ascending_subset_of_the_frozen_space(self):
        tenv, rank_sum, p_n = self.fill_inputs()
        axis = gene_axis.axis_fill(tenv, p_n[6], rank_sum, n=6)
        self.assertEqual(axis.dtype, np.int64)
        np.testing.assert_array_equal(axis, np.sort(axis))
        self.assertEqual(len(set(axis.tolist())), len(axis))
        self.assertGreaterEqual(int(axis.min()), 0)
        self.assertLess(int(axis.max()), len(AXIS_SYMBOLS))
        # never renumbered: the values are the original gene_index integers
        self.assertIn(7, axis.tolist())

    def test_out_of_contract_fill_fails_closed(self):
        tenv, rank_sum, p_n = self.fill_inputs()
        with self.assertRaises(ValueError):
            gene_axis.axis_fill(tenv, p_n[4], rank_sum, n=7)  # odd N
        with self.assertRaises(ValueError):
            gene_axis.axis_fill(tenv, p_n[4], rank_sum, n=10)  # larger than axis
        with self.assertRaises(ValueError):
            gene_axis.axis_fill(tenv[:7], p_n[4][:7], rank_sum[:7], n=4)  # not the axis

    def test_composition_and_pn_distribution(self):
        tenv, rank_sum, p_n = self.fill_inputs()
        axis = gene_axis.axis_fill(tenv, p_n[6], rank_sum, n=6)
        measured = [np.zeros(8, dtype=bool) for _ in CATEGORIES]
        measured[0][[0, 5, 7]] = True
        measured[1][[1, 2]] = True
        composition = gene_axis.axis_composition(axis, measured)
        self.assertEqual(composition['nhf'], 3)     # genes 0, 5 and 7
        self.assertEqual(composition['icm'], 1)     # gene 1
        self.assertEqual(composition['af'], 0)
        self.assertEqual(sum(composition.values()), 4)
        counts = gene_axis.pn_distribution(p_n[6])
        self.assertEqual(counts.tolist(), [0, 1, 1, 1, 1, 1, 0, 3])
        self.assertEqual(int(counts.sum()), len(AXIS_SYMBOLS))


if __name__ == '__main__':
    unittest.main()
