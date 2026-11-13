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
sys.path.insert(0, str(STEP / 'figure/script'))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
import build_gene_axis as driver
import render_gene_axis_rank_sum_histogram as renderer
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


def dataset_record(dataset_id, tier, isogenic=False, n_measured=8):
    record = {'dataset_id': dataset_id, 'tier': tier, 'n_measured': n_measured}
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
    'af': [('af_a', 1, False, IDENTITY)],
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


def write_fixture(root, *, drop_method=False, drop_derived=False):
    """Write a hermetic mini repository: axis, tenv ranking and six catalogs."""
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
            self.assertEqual(len(reader.composition()), 2 * (len(CATEGORIES) + 1 + 8 + 8))
            for text in manifest_strings(manifest):
                self.assertFalse(text.startswith('/'), f'absolute path in manifest: {text}')

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


class FigureTests(unittest.TestCase):
    """The evidence PNG is a function of its plotted-values CSV alone."""

    @staticmethod
    def synthetic_rows():
        rows = [{'record': 'figure', 'title': 'fixture'},
                {'record': 'panel', 'panel': 'a', 'panel_title': 'rank-sum histogram'}]
        for index in range(4):
            rows.append({'record': 'hist', 'bin_left': index, 'bin_right': index + 1,
                         'bin_center': index + 0.5, 'count': 10 + index,
                         'share': (10 + index) / 46})
        for axis, low, high, median in (('2000', 0.5, 2.0, 1.0),
                                        ('5000', 0.9, 3.0, 1.8)):
            rows.append({'record': 'interval', 'axis': axis, 'interval_low': low,
                         'interval_high': high, 'median': median})
        rows.append({'record': 'marker', 'inflection': False, 'inflection_n': '',
                     'note': 'no inflection'})
        return rows

    def test_renderer_rejects_incomplete_csv(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'bad.csv'
            path.write_text('record,panel,panel_title,title\nfigure,,,,x\n')
            with self.assertRaises(ValueError):
                renderer.load(path)

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


if __name__ == '__main__':
    unittest.main()
