"""Tests for the two-level rank-sum fusion and the released axes.

Everything is a hand-built 8-gene fixture except the real-file tests: the rule
document and the two committed figure PNGs are regenerated from their tracked
inputs and compared byte-for-byte, and the document's input wiring, its
source-table/excluded-AF contracts and its sample-size footnote are asserted
against the committed file. No release payload (h5ad, ``source/`` tree or
gitignored intermediate) is read.
The fusion, the axis fill rule, the readers and the document tables are
exercised at the public API / product boundary; one seam only (ADR 0066 /
spec). Determinism is tested by building the same fixture twice into scratch
candidate roots and comparing every digest.
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
from bases import deg_rank, gene_axis
from bases.deg_rank import normalized_rank
from bases.disease_deg import catalog_digest
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


def dataset_record(dataset_id, tier, isogenic=False, n_measured=8,
                   informative=True):
    record = {'dataset_id': dataset_id, 'tier': tier, 'n_measured': n_measured,
              'informative': informative}
    if not informative:
        record['informativeness'] = {'rule': 'fixture', 'n_distinct_effect': 1,
                                     'reason': 'all effect scores are tied'}
    if isogenic:
        record['isogenic_exception'] = {'applies': True, 'why': 'fixture isogenic pair'}
    return record


class ReExportContractTests(unittest.TestCase):
    """#186 contract: gene_axis re-exports the DEG/rank functions, never copies them."""

    def test_deg_rank_contract_is_re_exported_by_identity(self):
        for name in ('SymbolResolver', 'complete_rank', 'deg_effect', 'per_gene_deg',
                     'rank_scores', 'ALIGNMENT_CATEGORIES', 'PADJ_FLOOR',
                     'PADJ_THRESHOLD', 'normalized_rank', 'AXIS_SIZE'):
            self.assertIs(getattr(gene_axis, name), getattr(deg_rank, name),
                          f'bases.gene_axis.{name} is not bases.deg_rank.{name}')


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
        # Every tier-2 record is retained by the isogenic exception, so no tier
        # is excluded; the manifest must not call tier 2 excluded (ADR 0068).
        self.assertEqual(selected['excluded_tiers'], [])

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
        # unmeasured rank 7 / n_measured 6 = 7/6 (ADR 0066), never 1.0
        np.testing.assert_allclose(result['sum'], [1/6, 2/6, 3/6, 4/6, 5/6, 1.0, 7/6, 7/6])
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
    """Level 2: rank_sum over the six sources, ascending-is-better, pN, coverage."""

    @staticmethod
    def identity():
        return np.arange(1, len(AXIS_SYMBOLS) + 1, dtype=np.int64)

    def test_sum_of_six_normalized_ranks_and_gene_index_tie_break(self):
        ranks = [self.identity().copy() for _ in CATEGORIES]
        # one source puts gene 0 at rank 8 and gene 1 at rank 2; with the five
        # identity categories plus tenv-HVG the integer sums are 13 and 12.
        ranks[0] = np.array([8, 2, 3, 4, 5, 6, 7, 1], dtype=np.int64)
        measured = [np.ones(8, dtype=bool) for _ in CATEGORIES]
        result = gene_axis.fusion_rank_sum(
            ranks, measured, self.identity(), sizes=(2, 4))
        np.testing.assert_allclose(
            result['rank_sum'],
            np.array([13, 12, 18, 24, 30, 36, 42, 41]) / 8.0)
        # equal rank sums are ordered by gene_index in the fusion (fill-rule test)
        self.assertLess(result['rank_sum'][1], result['rank_sum'][0])

    def test_pn_boundaries_zero_and_six(self):
        ranks = [self.identity().copy() for _ in CATEGORIES]
        ranks[0][:2] = [2, 1]
        measured = [np.ones(8, dtype=bool) for _ in CATEGORIES]
        result = gene_axis.fusion_rank_sum(
            ranks, measured, self.identity(), sizes=(1, 2, 4))
        self.assertEqual(result['pn'][1].tolist(), [5, 1, 0, 0, 0, 0, 0, 0])
        self.assertEqual(result['pn'][2].tolist(), [6, 6, 0, 0, 0, 0, 0, 0])
        self.assertEqual(result['pn'][4].tolist(), [6, 6, 6, 6, 0, 0, 0, 0])
        self.assertEqual(int(result['pn'][2].max()), 6)
        self.assertEqual(int(result['pn'][2].min()), 0)

    def test_excluded_category_shrinks_the_pn_bound(self):
        # ADR 0068: with one of the five categories non-informative, four
        # categories + tenv-HVG make pN bounded by 5.
        ranks = [self.identity().copy() for _ in CATEGORIES[:4]]
        measured = [np.ones(8, dtype=bool) for _ in CATEGORIES[:4]]
        result = gene_axis.fusion_rank_sum(
            ranks, measured, self.identity(), categories=CATEGORIES[:4],
            sizes=(2,))
        self.assertEqual(result['n_sources'], 5)
        self.assertEqual(int(result['pn'][2].max()), 5)
        self.assertEqual(int(result['n_sources_measured'].max()), 5)
        # One rank per active category is required: four names, four ranks.
        with self.assertRaises(ValueError):
            gene_axis.fusion_rank_sum(
                ranks + [self.identity().copy()], measured + [np.ones(8, bool)],
                self.identity(), categories=CATEGORIES[:4], sizes=(2,))

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
        self.assertEqual(int(table.n_sources_measured[3]), 5)
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
        p_n = {4: np.array([0, 0, 0, 4, 0, 5, 6, 6]),
               6: np.array([4, 3, 2, 1, 5, 6, 6, 6])}
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
        p_n = np.array([6, 6, 6, 6, 0, 6, 6, 6])
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
        self.assertEqual(composition['nhf'], 3)     # genes 0, 5 and 7 are in G_6
        self.assertEqual(composition['icm'], 1)     # gene 1
        self.assertEqual(sum(composition.values()), 4)
        counts = gene_axis.pn_distribution(p_n[6])
        self.assertEqual(counts.tolist(), [0, 1, 1, 1, 1, 1, 3])
        self.assertEqual(int(counts.sum()), len(AXIS_SYMBOLS))


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

IDENTITY = np.arange(1, len(AXIS_SYMBOLS) + 1, dtype=np.int64)
SECOND = np.array([8, 1, 2, 3, 4, 5, 6, 7], dtype=np.int64)
FIXTURE_CATEGORIES = {
    'nhf': [('nhf_a', 1, False, IDENTITY)],
    'icm': [('icm_a', 1, False, IDENTITY)],
    'genetic_cardiomyopathy': [
        ('gcm_tier1', 1, False, IDENTITY),
        ('gcm_isogenic', 2, True, SECOND),
        ('gcm_dropped', 2, False, SECOND),
    ],
    'cardiotox_TOP2': [('top2_a', 2, False, IDENTITY)],
    'cardiotox_other': [('other_a', 2, False, IDENTITY)],
}
FIXTURE_AXIS_SIZES = (4, 6)
FIXTURE_TENV_RANK = np.array([8, 7, 6, 5, 4, 3, 2, 1], dtype=np.int64)


def write_fixture(root, *, drop_method=False, drop_derived=False,
                  noninformative=()):
    """Write a hermetic mini repository: axis, tenv ranking and five catalogs."""
    root = Path(root)
    step = root / 'stage2/01_bases'
    genes = axis_fixture()
    genes_path = step / 'output/gene-identity/genes.csv'
    genes_path.parent.mkdir(parents=True, exist_ok=True)
    genes.to_csv(genes_path, index=False, lineterminator='\n')
    tenv = pd.DataFrame({
        'gene_index': np.arange(len(AXIS_SYMBOLS)),
        'log_mean': np.zeros(len(AXIS_SYMBOLS)),
        'log_dispersion': np.zeros(len(AXIS_SYMBOLS)),
        'normalized_dispersion': -FIXTURE_TENV_RANK.astype(float),
        'rank': FIXTURE_TENV_RANK,
    })
    tenv_path = step / 'output/tenv-hvg/descriptive.csv'
    tenv_path.parent.mkdir(parents=True, exist_ok=True)
    tenv.to_csv(tenv_path, index=False, line_terminator='\n') if False else tenv.to_csv(
        tenv_path, index=False, lineterminator='\n')
    disease = step / 'output/disease-deg'
    for category, records in FIXTURE_CATEGORIES.items():
        directory = disease / category
        directory.mkdir(parents=True, exist_ok=True)
        canon = category not in ('cardiotox_TOP2', 'cardiotox_other')
        catalog_records, outputs = [], {}
        for dataset_id, tier, isogenic, rank in records:
            informative = dataset_id not in noninformative
            rank_name = f'{dataset_id}.rank.csv'
            frame = pd.DataFrame({'gene_index': np.arange(len(AXIS_SYMBOLS)),
                                  'rank': rank,
                                  'rank_norm': normalized_rank(rank, len(AXIS_SYMBOLS))})
            frame.to_csv(directory / rank_name, index=False, lineterminator='\n')
            outputs[rank_name] = {'sha256': hashlib.sha256(
                (directory / rank_name).read_bytes()).hexdigest()}
            record = {
                'dataset_id': dataset_id, 'tier': tier,
                'evidence_type': 'cm_pseudobulk' if canon else 'bulk',
                'n_measured': len(AXIS_SYMBOLS),
                'n_significant': 2, 'n_case': 3, 'n_control': 3,
                'informative': informative,
                'informativeness': {'rule': 'fixture', 'n_distinct_effect': 2,
                                    'reason': None},
                'deg': f'{dataset_id}.deg.csv', 'rank': rank_name,
                'method': {'model': 'fixture'},
                'inputs': {'fixture_input.csv': {
                    'sha256': hashlib.sha256(b'fixture-input').hexdigest()}},
            }
            if isogenic:
                record['isogenic_exception'] = {'applies': True}
            if canon:
                record['derived_inputs'] = [{
                    'relative': 'local/deg/fixture.csv',
                    'sha256': hashlib.sha256(b'fixture-derived').hexdigest()}]
            if drop_method and dataset_id == 'nhf_a':
                record.pop('method')
            if drop_derived and dataset_id == 'nhf_a':
                record.pop('derived_inputs')
            catalog_records.append(record)
        catalog = {
            'schema': 'stage2-disease-deg-catalog-v1',
            'category': {'nhf': 'nHF'}.get(category, category),
            'directory': category, 'axis_size': len(AXIS_SYMBOLS),
            'tier': 1 if canon else 2,
            'evidence_type': 'cm_pseudobulk' if canon else 'bulk',
            'datasets': catalog_records,
        }
        if canon:
            catalog['rank_contract'] = 'stage2-disease-deg-v1'
        catalog['catalog_sha256'] = catalog_digest(catalog)
        (directory / 'catalog.json').write_text(
            json.dumps(catalog, indent=2, sort_keys=True) + '\n')
        (directory / 'manifest.json').write_text(json.dumps({
            'schema': 'stage2-disease-deg-release-v1',
            'catalog_sha256': catalog['catalog_sha256'],
            'outputs': dict(sorted(outputs.items())),
        }, indent=2, sort_keys=True) + '\n')
    return driver.default_config(
        root=root,
        output=step / 'output/gene-axis',
        figure_data=step / 'figure/data/gene_axis_rank_sum_histogram.csv',
        axis_sizes=FIXTURE_AXIS_SIZES)


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


class ReleaseBuildTests(unittest.TestCase):
    """Candidate-then-replace publication, replay determinism and the reader seam."""

    def test_publish_products_and_reader_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            result = driver.publish(config)
            self.assertEqual(result['n_candidates'], len(AXIS_SYMBOLS))
            reader = GeneAxis(config.output)
            manifest = reader.manifest()
            self.assertEqual(manifest['schema'], gene_axis.SCHEMA)
            self.assertEqual(reader.axis_sizes(), FIXTURE_AXIS_SIZES)
            table = reader.candidate_table()
            self.assertEqual(list(table.columns),
                             list(gene_axis.candidate_columns(FIXTURE_AXIS_SIZES)))
            self.assertEqual(len(table), len(AXIS_SYMBOLS))
            np.testing.assert_array_equal(reader.axis(4), np.array([0, 1, 6, 7]))
            np.testing.assert_array_equal(reader.axis(6), np.array([2, 3, 4, 5, 6, 7]))
            self.assertFalse(set(reader.axis(4)) <= set(reader.axis(6)))
            self.assertEqual(manifest['intersection']['n_intersection'], 2)
            self.assertEqual(manifest['intersection']['sizes'], [4, 6])
            self.assertEqual(manifest['categories']['genetic_cardiomyopathy']
                             ['winning_tier'], 1)
            self.assertEqual(
                manifest['categories']['genetic_cardiomyopathy']['retained_by_exception'],
                ['gcm_isogenic'])
            self.assertNotIn('gcm_dropped',
                             [record['dataset_id'] for record in
                              manifest['categories']['genetic_cardiomyopathy']['datasets']])
            # 2 axes x (5 categories + tenv half + 7 pN values + 7 candidate pN)
            self.assertEqual(len(reader.composition()), 2 * (len(CATEGORIES) + 1 + 7 + 7))
            coverage = reader.coverage()
            self.assertEqual(list(coverage.columns), list(gene_axis.COVERAGE_COLUMNS))
            self.assertEqual(coverage.N.tolist(), list(FIXTURE_AXIS_SIZES))
            self.assertEqual(coverage['%tenv_hvg_gN'].tolist(), [50.0, 100.0])
            self.assertIn('source_topn_coverage.csv',
                          ' '.join(manifest['outputs'].keys()))
            self.assertIn('gene_axis_source_topn_coverage.csv',
                          ' '.join(manifest['outputs'].keys()))
            coverage_figure = pd.read_csv(config.source_coverage_figure)
            self.assertEqual(list(coverage_figure.columns),
                             list(driver.SOURCE_COVERAGE_FIGURE_COLUMNS))
            self.assertEqual(len(coverage_figure),
                             2 + len(FIXTURE_AXIS_SIZES) * (len(CATEGORIES) + 1))
            bars = coverage_figure[coverage_figure.record == 'bar']
            self.assertTrue((bars['count'] % 1 == 0).all())
            self.assertTrue((bars['count'] >= 0).all())
            self.assertTrue(bars['rate'].between(0.0, 100.0).all())
            self.assertEqual(bars[bars.source == 'tenv_hvg'].set_index('series')
                             ['count'].to_dict(), {'G_4': 2, 'G_6': 6})
            panels = coverage_figure[coverage_figure.record == 'panel']
            self.assertEqual(panels['panel_title'].tolist(),
                             ['Covered Gene Count', 'Top-N Coverage Rate'])
            self.assertEqual(manifest['coverage']['figure']['panels'],
                             [dict(panel) for panel in
                              driver.SOURCE_COVERAGE_FIGURE_PANELS])
            for text in manifest_strings(manifest):
                self.assertFalse(text.startswith('/'), f'absolute path in manifest: {text}')

    def test_noninformative_category_is_excluded_from_fusion_and_pn(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary), noninformative=('icm_a',))
            driver.publish(config)
            manifest = GeneAxis(config.output).manifest()
            self.assertEqual(manifest['fusion']['n_sources'], 5)
            self.assertEqual(manifest['fusion']['excluded_categories'], ['icm'])
            self.assertEqual(manifest['fusion']['sources'],
                             [category for category in CATEGORIES if category != 'icm']
                             + ['tenv_hvg'])
            self.assertFalse(manifest['categories']['icm']['informative'])
            self.assertTrue(manifest['categories']['icm']['excluded_from_fusion'])
            self.assertEqual(manifest['categories']['icm']['excluded_datasets'][0]
                             ['dataset_id'], 'icm_a')
            self.assertEqual(manifest['release']['pN_bounds'], [0, 5])
            table = GeneAxis(config.output).candidate_table()
            self.assertTrue(table['icm_rank'].isna().all())
            for size in FIXTURE_AXIS_SIZES:
                self.assertLessEqual(int(table[f'p{size}'].max()), 5)
                self.assertLessEqual(
                    int(manifest['axes'][str(size)]['pN_distribution'].get('6', 0)), 0)
            composition = GeneAxis(config.output).composition()
            icm = composition[(composition.record == 'composition')
                              & (composition.group == 'icm')]
            self.assertTrue((icm.n_genes == 0).all())
            self.assertTrue(icm.note.str.contains('excluded').all())
            coverage = GeneAxis(config.output).coverage()
            self.assertTrue(coverage['%icm_gN'].isna().all())
            # Every other category still fuses: four + tenv.
            self.assertEqual(driver.verify(config)['result'], 'PASS')

    def test_publish_preserves_the_retired_gene_axis_file_set(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary))
            config.output.mkdir(parents=True, exist_ok=True)
            legacy = config.output / 'fused_rank.csv'
            legacy.write_text('retired product\n')
            driver.publish(config)
            self.assertTrue(legacy.is_file())  # removal is #192's job, not #190's

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

    def test_uncontracted_catalog_shape_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary), drop_method=True)
            with self.assertRaises(ValueError):
                driver.publish(config)
        with tempfile.TemporaryDirectory() as temporary:
            config = write_fixture(Path(temporary), drop_derived=True)
            with self.assertRaises(ValueError):
                driver.publish(config)


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class CoverageTests(unittest.TestCase):
    """Each source's top-N coverage is a pure function of the axes and its rank."""

    @staticmethod
    def ranks():
        return {'tenv_hvg': np.array([8, 7, 6, 5, 4, 3, 2, 1]),
                'nhf': np.array([1, 2, 3, 4, 5, 6, 7, 8]),
                'icm': None,
                'genetic_cardiomyopathy': np.array([1, 2, 3, 4, 5, 6, 7, 8]),
                'cardiotox_TOP2': np.array([2, 1, 4, 3, 6, 5, 8, 7]),
                'cardiotox_other': np.array([4, 3, 2, 1, 8, 7, 6, 5])}

    def test_percentages_and_empty_excluded_column(self):
        axes = {4: np.array([0, 1, 6, 7]), 6: np.array([2, 3, 4, 5, 6, 7])}
        table = gene_axis.source_topn_coverage(axes, self.ranks(), sizes=(4, 6))
        self.assertEqual(list(table.columns), list(gene_axis.COVERAGE_COLUMNS))
        self.assertEqual(table.N.tolist(), [4, 6])
        # tenv top-4 = genes 4..7 and G_4 holds 6,7 -> 2/4 = 50 %.
        self.assertEqual(table['%tenv_hvg_gN'][0], 50.0)
        # tenv top-6 = genes 2..7, exactly G_6 -> 100 %.
        self.assertEqual(table['%tenv_hvg_gN'][1], 100.0)
        # nhf ranks identity: top-4 = genes 0..3, G_4 holds 0,1 -> 50 %.
        self.assertEqual(table['%nhf_gN'][0], 50.0)
        self.assertTrue(table['%icm_gN'].isna().all())

    def test_covered_counts_are_the_percentage_numerators(self):
        axes = {4: np.array([0, 1, 6, 7]), 6: np.array([2, 3, 4, 5, 6, 7])}
        ranks = self.ranks()
        counts = gene_axis.source_topn_covered_counts(axes, ranks, sizes=(4, 6))
        table = gene_axis.source_topn_coverage(axes, ranks, sizes=(4, 6))
        self.assertEqual(list(counts.columns), list(gene_axis.COUNT_COLUMNS))
        self.assertEqual(counts['tenv_hvg'].tolist(), [2, 6])
        self.assertEqual(counts['nhf'].tolist(), [2, 4])
        self.assertTrue(counts['icm'].isna().all())
        for index, size in enumerate((4, 6)):
            for name in ('tenv_hvg', 'nhf', 'genetic_cardiomyopathy',
                         'cardiotox_TOP2', 'cardiotox_other'):
                count = int(counts[name].iloc[index])
                self.assertGreaterEqual(count, 0)
                self.assertLessEqual(count, size)
                self.assertEqual(table[f'%{name}_gN'].iloc[index],
                                 round(100.0 * count / size, 2))

    def test_out_of_contract_inputs_fail_closed(self):
        axes = {4: np.array([0, 1, 6, 7])}
        bad = self.ranks()
        bad['tenv_hvg'] = np.array([1, 1, 2, 3, 4, 5, 6, 7])
        with self.assertRaises(ValueError):
            gene_axis.source_topn_coverage(axes, bad, sizes=(4,))
        with self.assertRaises(ValueError):
            gene_axis.source_topn_coverage(axes, self.ranks(), sizes=(5,))
        with self.assertRaises(ValueError):
            gene_axis.source_topn_coverage(
                axes, {'unknown': np.arange(1, 9)}, sizes=(4,))


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
    for category in gene_axis.CATEGORIES:
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
                for category in gene_axis.CATEGORIES}

    def test_committed_docx_replays_from_its_inputs(self):
        coverage = STEP / 'output/gene-axis/source_topn_coverage.csv'
        figure = STEP / 'figure/gene_axis_rank_sum_histogram.png'
        catalogs = self.catalogs()
        committed = STEP / 'output/gene-axis/gene_axis_fill_rule.docx'
        for path in (coverage, figure, committed, *catalogs.values()):
            self.assertTrue(path.is_file(), f'missing docx input or output: {path}')
        import gene_axis_source_metadata as metadata
        import make_gene_axis_docx as docx_builder
        with tempfile.TemporaryDirectory() as temporary:
            rendered = Path(temporary) / 'gene_axis_fill_rule.docx'
            docx_builder.build(coverage, figure, rendered,
                               catalogs=catalogs, metadata=metadata)
            self.assertEqual(sha256_file(rendered), sha256_file(committed))

    def test_build_reads_the_five_catalog_paths_it_is_given(self):
        """A changed catalog ``n_significant`` must reach the rendered cell."""
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
                    dataset['n_significant'] = 12345
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
            catalog['datasets'][0]['n_significant'] = 424242
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
        coverage = STEP / 'output/gene-axis/source_topn_coverage.csv'
        import make_gene_axis_docx as docx_builder
        table = docx_builder.read_coverage(coverage)
        self.assertEqual(list(table.columns), list(docx_builder.COVERAGE_COLUMNS))

    def test_source_table_rows_match_catalogs(self):
        """The first table is the released source dataset table, row for row."""
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
        self.assertEqual(len(table.rows), 1 + 19, 'expected 19 source rows + header')
        expected = []
        for category in gene_axis.CATEGORIES:
            catalog = json.loads(self.catalogs()[category].read_text(encoding='utf-8'))
            expected.extend((catalog['category'], dataset) for dataset in catalog['datasets'])
        self.assertEqual(len(expected), 18)
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
                self.assertEqual(int(cells[6].replace(',', '')), dataset['n_significant'],
                                 'CM DEG count column')
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

    def test_af_excluded_table_contract(self):
        """The second table is the audited-and-excluded AF sources block."""
        from docx import Document
        import make_gene_axis_docx as docx_builder
        committed = STEP / 'output/gene-axis/gene_axis_fill_rule.docx'
        self.assertTrue(committed.is_file(), f'missing committed docx: {committed}')
        document = Document(str(committed))
        self.assertGreaterEqual(len(document.tables), 2, 'the document has no AF table')
        table = document.tables[1]
        header = [cell.text for cell in table.rows[0].cells]
        self.assertEqual(header, list(docx_builder.AF_TABLE_COLUMNS))
        self.assertEqual(len(table.rows), 1 + 5, 'expected 5 excluded-AF rows + header')
        rows = [[cell.text for cell in row.cells] for row in table.rows[1:]]
        self.assertEqual([row[0] for row in rows],
                         ['GSE255612', 'GSE138252', 'GSE238242', 'GSE271229',
                          'dbGaP phs003644'])
        expected_deg = {
            'GSE255612': ('0 in this pipeline', '0.7359', '2 in the paper',
                          'LINC01257', 'ATRNL1'),
            'GSE138252': ('2 (limma-voom)', '6 (edgeR QL)'),
            'GSE238242': ('0 (min padj 0.217)',),
            'GSE271229': ('not testable',),
            'dbGaP phs003644': ('not downloadable',),
        }
        expected_reason = {
            'GSE255612': 'no AF-specific signal',
            'GSE138252': 'confound',
            'GSE238242': 'no signal at n=4/3',
            'GSE271229': 'no cell-type annotation',
            'dbGaP phs003644': 'controlled access',
        }
        for row in rows:
            with self.subTest(source=row[0]):
                for fragment in expected_deg[row[0]]:
                    self.assertIn(fragment, row[2], 'actual CM DEG count cell')
                self.assertIn(expected_reason[row[0]], row[3], 'exclusion reason cell')


if __name__ == '__main__':
    unittest.main()
