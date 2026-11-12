"""Fixture-only tests for the gene-axis evidence definitions and its seams.

No real data is read: the axis, expression, DEG, ortholog and plotted-value
tables below are hand-built fixtures. The only repository files opened are the
tracked evidence artifacts (the committed figure CSV/PNG and the committed
manifest) which the replay and leakage tests must verify without touching
``output/`` products.
"""
import builtins
import csv
import hashlib
import io
import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

EXECUTION = Path(__file__).resolve().parents[1]
STEP = EXECUTION.parent
sys.path.insert(0, str(EXECUTION))
sys.path.insert(0, str(STEP / 'figure/script'))
sys.path.insert(0, str(STEP / 'output/api'))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
import build_gene_axis as driver
import render_gene_axis_source_summary as renderer
from bases.gene_axis import (
    AXIS_SIZE, SymbolResolver, align_to_gene_index, collapse_expression, complete_rank,
    deg_count_agreement, deg_effect, detect_inflection, direction_composition, fuse_ranks,
    fused_rank_stats, hvg_membership, kneedle_knee, normalized_dispersion, normalized_rank,
    overlap_metrics, per_gene_deg, rank_scores, second_difference_peak, tenv_non_hvg_curve,
    transcriptome_hvg_sets,
)

AXIS_SYMBOLS = ['G0', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7']


def axis_fixture():
    return pd.DataFrame({
        'gene_index': np.arange(len(AXIS_SYMBOLS)),
        'symbol': AXIS_SYMBOLS,
        'hgnc_id': [f'HGNC:{index + 1}' for index in range(len(AXIS_SYMBOLS))],
        'entrez_id': [str(index + 1) for index in range(len(AXIS_SYMBOLS))],
        'ensembl_id': [f'ENSG0000000000{index + 1}' for index in range(len(AXIS_SYMBOLS))],
    })


def crosswalk_fixture():
    return pd.DataFrame([
        {'dataset_id': 'dtoxs', 'source_gene_id': 'OLD_G0', 'symbol': 'G0',
         'hgnc_id': 'HGNC:1', 'match_type': 'previous'},
        {'dataset_id': 'dtoxs', 'source_gene_id': 'ALIAS_G1', 'symbol': 'G1',
         'hgnc_id': 'HGNC:2', 'match_type': 'alias'},
        {'dataset_id': 'dtoxs', 'source_gene_id': 'AMB', 'symbol': 'X',
         'hgnc_id': 'HGNC:3', 'match_type': 'alias'},
        {'dataset_id': 'dtoxs', 'source_gene_id': 'AMB', 'symbol': 'Y',
         'hgnc_id': 'HGNC:4', 'match_type': 'alias'},
        {'dataset_id': 'dtoxs', 'source_gene_id': 'OFF', 'symbol': 'Z',
         'hgnc_id': 'HGNC:999', 'match_type': 'previous'},
        {'dataset_id': 'dtoxs', 'source_gene_id': 'ENSG00000000005', 'symbol': 'G4',
         'hgnc_id': 'HGNC:5', 'match_type': 'ensembl'},
    ])


def hgnc_fixture():
    return pd.DataFrame([
        {'hgnc_id': 'HGNC:1', 'symbol': 'G0', 'alias_symbol': 'HG0|HG0B',
         'prev_symbol': '', 'ensembl_gene_id': 'ENSG00000000001'},
        {'hgnc_id': 'HGNC:6', 'symbol': 'G5', 'alias_symbol': 'HG5',
         'prev_symbol': 'OLD5', 'ensembl_gene_id': 'ENSG00000000006'},
        {'hgnc_id': 'HGNC:999', 'symbol': 'OFFAXIS', 'alias_symbol': '',
         'prev_symbol': '', 'ensembl_gene_id': 'ENSG000000009999'},
    ])


class AlignmentTests(unittest.TestCase):
    def setUp(self):
        self.resolver = SymbolResolver(axis_fixture(), crosswalk_fixture().to_dict('records'),
                                       hgnc_fixture().to_dict('records'))

    def test_precedence_and_categories(self):
        cases = {
            'G0': (0, 'canonical'),
            'OLD_G0': (0, 'crosswalk_alias'),
            'ALIAS_G1': (1, 'crosswalk_alias'),
            'HG5': (5, 'hgnc_alias'),
            'OLD5': (5, 'hgnc_alias'),
            'ENSG00000000005': (4, 'ensembl'),
            'ENSG00000000005.7': (4, 'ensembl'),
            'AMB': (-1, 'ambiguous'),
            'OFF': (-1, 'offaxis'),
            'UNKNOWN_SYMBOL': (-1, 'unmapped'),
            'ENSG000000009999': (-1, 'offaxis'),
        }
        for symbol, expected in cases.items():
            self.assertEqual(self.resolver.resolve(symbol), expected, symbol)

    def test_loss_summary_counts_every_symbol(self):
        symbols = ['G0', 'G1', 'AMB', 'OFF', 'UNKNOWN_SYMBOL', 'HG5']
        indices, categories = self.resolver.resolve_many(symbols)
        summary = self.resolver.loss_summary(indices, categories)
        self.assertEqual(summary['n_source_symbols'], len(symbols))
        self.assertEqual(summary['canonical'], 2)
        self.assertEqual(summary['hgnc_alias'], 1)
        self.assertEqual(summary['ambiguous'], 1)
        self.assertEqual(summary['offaxis'], 1)
        self.assertEqual(summary['unmapped'], 1)
        self.assertEqual(summary['n_resolved'], 3)
        self.assertEqual(summary['n_axis_genes_resolved'], 3)


class RankTests(unittest.TestCase):
    def test_effect_score_and_na_padj(self):
        effect = deg_effect(np.array([2.0, -2.0, 1.0, 3.0]), np.array([1e-4, 0.5, np.nan, 0.0]))
        self.assertAlmostEqual(effect[0], 2.0 * 4.0)
        self.assertAlmostEqual(effect[1], 2.0 * -np.log10(0.5))
        self.assertEqual(effect[2], 0.0)   # NA padj is non-significant, not dropped
        self.assertAlmostEqual(effect[3], 3.0 * 300.0)  # floored at 1e-300

    def test_multi_contrast_collapse_by_max(self):
        per_gene = per_gene_deg(
            gene_index=np.array([0, 0, 1, 1]),
            log2fc=np.array([0.5, -2.0, 1.0, 0.2]),
            padj=np.array([0.04, 0.001, 0.5, np.nan]),
            contrast=np.array(['a', 'b', 'a', 'b']))
        winning = per_gene.set_index('gene_index')
        self.assertEqual(winning.loc[0, 'contrast'], 'b')
        self.assertAlmostEqual(winning.loc[0, 'effect'], 2.0 * 3.0)
        self.assertEqual(winning.loc[0, 'direction'], -1)
        self.assertEqual(winning.loc[0, 'n_contrasts_supporting'], 2)
        self.assertEqual(winning.loc[1, 'contrast'], 'a')
        self.assertEqual(winning.loc[1, 'direction'], 0)   # winning padj is not significant
        self.assertEqual(winning.loc[1, 'n_contrasts_supporting'], 0)

    def test_unmeasured_gets_worst_rank_and_one(self):
        measured = np.array([0, 2, 4])
        effect = np.array([1.0, 3.0, 2.0])
        rank = complete_rank(measured, rank_scores(measured, effect), axis_size=8)
        self.assertEqual(rank.tolist(), [3, 4, 1, 4, 2, 4, 4, 4])
        norm = normalized_rank(rank, 3)
        self.assertEqual(norm[1], 1.0)                      # unmeasured
        self.assertEqual(norm[0], 1.0)
        self.assertEqual(norm[2], 1 / 3)

    def test_rank_ties_break_on_gene_index(self):
        rank = rank_scores(np.array([4, 2, 9, 1]), np.array([1.0, 1.0, 1.0, 1.0]))
        self.assertEqual(rank.tolist(), [3, 2, 4, 1])

    def test_hvg_recipe_ranks_variable_genes_first(self):
        values = np.array([
            [1.0, 1.0, 1.0, 1.0],
            [0.0, 0.0, 10.0, 10.0],
            [2.0, 2.1, 1.9, 2.0],
            [1.0, 3.0, 5.0, 7.0],
        ])
        table = normalized_dispersion(values, np.arange(4))
        self.assertEqual(int(table.sort_values('rank').gene_index.iloc[0]), 1)
        self.assertEqual(sorted(table['rank'].tolist()), [1, 2, 3, 4])


class FusionTests(unittest.TestCase):
    def test_unweighted_fusion_tie_break_and_animal_sensitivity(self):
        n = 4
        ranks = {
            'tenv': np.array([1, 2, 3, 4]),
            'human': np.array([1, 2, 3, 4]),
            'gse160516': np.array([1, 2, 3, 4]),
            'gse108940': np.array([2, 1, 3, 4]),
        }
        measured = {name: 4 for name in ranks}
        fused = fuse_ranks(ranks, measured, axis_size=n)
        self.assertEqual(fused.rank_sum.tolist(), [1.25, 1.75, 3.0, 4.0])
        self.assertEqual(fused.fused_rank.tolist(), [1, 2, 3, 4])
        # animal sources at 0.5: 0.25 + 0.25 + 0.5*0.25 + 0.5*0.5 = 0.875
        self.assertAlmostEqual(fused.rank_sum_animal_half.iloc[0], 0.875)
        self.assertAlmostEqual(fused.rank_sum_animal_half.iloc[1], 1.375)
        self.assertEqual(fused.n_sources_measured.tolist(), [4, 4, 4, 4])

    def test_unmeasured_counts_and_tie_break(self):
        ranks = {
            'a': np.array([1, 2, 3]),          # 2 measured
            'b': np.array([1, 1, 3]),          # 2 measured
        }
        fused = fuse_ranks(ranks, {'a': 2, 'b': 2}, axis_size=3)
        self.assertEqual(fused.rank_sum.tolist(), [1.0, 1.5, 2.0])
        self.assertEqual(fused.n_sources_measured.tolist(), [2, 2, 0])
        self.assertEqual(fused.fused_rank.tolist(), [1, 2, 3])  # tie on gene_index
        self.assertEqual(fused.n_sources_deg_significant.tolist(), [0, 0, 0])

    def test_fusion_is_deterministic_and_order_independent(self):
        rng = np.random.default_rng(17)
        ranks = {name: rng.permutation(6) + 1 for name in ('a', 'b', 'c')}
        measured = {'a': 4, 'b': 5, 'c': 3}
        first = fuse_ranks(ranks, measured, axis_size=6)
        second = fuse_ranks(dict(reversed(list(ranks.items()))), measured, axis_size=6)
        columns = sorted(first.columns)
        pd.testing.assert_frame_equal(first[columns], second[columns])
        pd.testing.assert_frame_equal(first, fuse_ranks(ranks, measured, axis_size=6))

    def test_overlap_metrics(self):
        fused_rank = np.array([3, 1, 2])   # gene-indexed: gene 1 is fused rank 1
        metrics = overlap_metrics({0, 1}, {1, 2}, fused_rank)
        self.assertEqual(metrics['n_intersection'], 1)
        self.assertAlmostEqual(metrics['jaccard'], 1 / 3)
        self.assertEqual(metrics['fused_rank_min'], 1)

    def test_tenv_non_hvg_curve(self):
        ranked_genes = np.array([2, 0, 1])   # fused order: genes 2, 0, 1
        tenv_rank = np.array([3, 1, 2])      # gene-indexed
        curve = tenv_non_hvg_curve(ranked_genes, tenv_rank, ks=(1, 2))
        self.assertEqual(curve.n.tolist(), [1, 2, 3])
        self.assertEqual(curve['n_non_hvg_k1'].tolist(), [1, 2, 2])
        self.assertEqual(curve['n_non_hvg_k2'].tolist(), [0, 1, 1])
        self.assertTrue(hvg_membership(tenv_rank, 1).tolist() == [False, True, False])


class IndexSpaceTests(unittest.TestCase):
    """The tracked fused table is sorted by fused_rank, never by gene_index.

    These fixtures make fused-rank order differ from gene_index order and pin
    that every consumer answers the fused-order question, not the row-position
    one (the bug class this round fixes).
    """

    def test_align_to_gene_index_scatters_sorted_rows(self):
        # Rows in fused-rank order hold gene indices [1, 3, 0, 4, 2].
        row_gene_index = np.array([1, 3, 0, 4, 2])
        row_fused_rank = np.array([1, 2, 3, 4, 5])
        scattered = align_to_gene_index(row_fused_rank, row_gene_index, axis_size=5)
        self.assertEqual(scattered.tolist(), [3, 1, 5, 2, 4])
        with self.assertRaises(ValueError):
            align_to_gene_index(np.array([1, 2]), np.array([0, 1]), axis_size=5)
        with self.assertRaises(ValueError):
            align_to_gene_index(np.array([1, 2]), np.array([0, 0]), axis_size=2)

    def test_tenv_curve_uses_fused_order_not_gene_index(self):
        tenv_rank = np.array([1, 2, 3, 4, 5])   # gene 0 is tenv's top HVG
        fused_rank = np.array([3, 1, 5, 2, 4])  # fused order genes [1, 3, 0, 4, 2]
        fused_order = np.argsort(fused_rank)
        curve = tenv_non_hvg_curve(fused_order, tenv_rank, ks=(1, 2))
        self.assertEqual(curve['n_non_hvg_k1'].tolist(), [1, 2, 2, 3, 4])
        self.assertEqual(curve['n_non_hvg_k2'].tolist(), [0, 1, 1, 2, 3])
        # The old row-position call (gene_index order) must give a different curve.
        wrong = tenv_non_hvg_curve(np.arange(5), tenv_rank, ks=(1,))
        self.assertEqual(wrong['n_non_hvg_k1'].tolist(), [0, 1, 2, 3, 4])

    def test_overlap_rank_stats_use_gene_indexed_fused_rank(self):
        fused_rank = np.array([3, 1, 2])
        own = overlap_metrics({0}, {0}, fused_rank)
        self.assertEqual(own['fused_rank_min'], 3)
        self.assertEqual(own['fused_rank_max'], 3)
        self.assertEqual(fused_rank_stats({1}, fused_rank)['fused_rank_min'], 1)
        self.assertEqual(fused_rank_stats({0, 2}, fused_rank)['fused_rank_median'], 2.5)
        # A row-position array would report 2 for gene 0, not its true rank 3.
        positional = np.array([1, 2, 3])
        self.assertEqual(overlap_metrics({0}, {0}, positional)['fused_rank_min'], 1)

    def test_direction_composition_rules_select_different_genes(self):
        fused_rank = np.array([1, 2, 3, 4, 5, 6])
        direction = np.array([1, 1, -1, -1, 1, -1])
        max_abs_effect = np.array([0.0, 0.0, 0.0, 0.0, 9.0, 8.0])
        frame = direction_composition(fused_rank, direction, max_abs_effect, sizes=(2, 6))
        symmetric = frame[frame.rule == 'symmetric'].set_index('n')
        weighted = frame[frame.rule == 'max_abs_log2fc'].set_index('n')
        self.assertEqual((symmetric.loc[2, 'n_up'], symmetric.loc[2, 'n_down']), (2, 0))
        self.assertEqual((weighted.loc[2, 'n_up'], weighted.loc[2, 'n_down']), (1, 1))
        self.assertAlmostEqual(symmetric.loc[2, 'up_share'], 1.0)
        self.assertAlmostEqual(weighted.loc[2, 'up_share'], 0.5)
        self.assertEqual((symmetric.loc[6, 'n_up'], weighted.loc[6, 'n_up']), (3, 3))

    def test_direction_composition_uses_fused_rank_values_not_row_order(self):
        fused_rank = np.array([3, 1, 5, 2, 4])
        direction = np.array([1, -1, 1, -1, 1])
        max_abs_effect = np.zeros(5)
        frame = direction_composition(fused_rank, direction, max_abs_effect, sizes=(2,))
        symmetric = frame[frame.rule == 'symmetric'].iloc[0]
        self.assertEqual((symmetric.n_up, symmetric.n_down), (0, 2))  # genes 1, 3
        wrong = direction_composition(np.arange(1, 6), direction, max_abs_effect, sizes=(2,))
        wrong_row = wrong[wrong.rule == 'symmetric'].iloc[0]
        self.assertEqual((wrong_row.n_up, wrong_row.n_down), (1, 1))  # row positions 0, 1

    def test_deg_count_tie_break_is_fused_rank(self):
        fused_rank = np.array([4, 5, 1, 2, 3, 6])
        count = np.array([1, 1, 0, 1, 0, 0])
        result = deg_count_agreement(fused_rank, count, sizes=(2,))
        self.assertEqual(result['top_n_agreement']['2']['n_overlap'], 1)
        self.assertAlmostEqual(result['top_n_agreement']['2']['share'], 0.5)
        # The frozen gene_index tie-break is a different, weaker ranking.
        gene_order = np.lexsort((np.arange(6), -count))
        fused_order = np.argsort(fused_rank)
        self.assertEqual(len(set(fused_order[:2]) & set(gene_order[:2])), 0)

    def test_transcriptome_hvg_sets_project_top_k_onto_axis(self):
        hvg_rank = np.array([1, 2, 3, 4])
        row_gene_index = np.array([-1, 0, 1, 2])   # the most variable row is off-axis
        sets = transcriptome_hvg_sets(hvg_rank, row_gene_index, ks=(1, 2, 4), axis_size=3)
        self.assertEqual(sets[1]['n_total'], 1)
        self.assertEqual(sets[1]['n_on_axis'], 0)
        self.assertFalse(sets[1]['set'].any())
        self.assertEqual(sets[2]['n_on_axis'], 1)
        self.assertEqual(sets[4]['n_total'], 4)
        self.assertEqual(sets[4]['n_on_axis'], 3)
        self.assertAlmostEqual(sets[2]['n_on_axis'] / sets[2]['n_total'], 0.5)

    def test_source_summary_hvg_share_and_density(self):
        spec = {'id': 'x', 'accession': 'X', 'category': 'disease_data',
                'species': 'Homo sapiens', 'tissue_model': 't', 'contrast': 'c',
                'method': 'm', 'ortholog_chain': 'none', 'fused': True,
                'exploratory': False, 'supplementary_only': False, 'notes': ''}
        axis_set = np.zeros(8, dtype=bool)
        axis_set[[0, 1, 2, 3, 4]] = True
        result = {'measured': 6, 'n_axis_expression': 7, 'n_samples': 2,
                  'n_contrasts': 1, 'per_gene': None, 'metadata': {},
                  'hvg_sets': {2000: {'set': axis_set, 'n_total': 8, 'n_on_axis': 5}}}
        row = driver._source_summary_row(spec, result, axis_size=8)
        self.assertEqual(row['n_hvg'], 8)
        self.assertEqual(row['n_hvg_total'], 8)
        self.assertEqual(row['n_hvg_on_axis'], 5)
        self.assertAlmostEqual(row['hvg_share'], 5 / 8)
        self.assertAlmostEqual(row['hvg_density'], 5 / 7)


class InflectionTests(unittest.TestCase):
    @staticmethod
    def knee_curve(knee=750, length=2000):
        x = np.arange(1, length + 1, dtype=float)
        y = np.minimum(x * 0.001, 0.6 + 0.0002 * (x - 1))
        return x, y

    def test_synthetic_knee_is_detected(self):
        x, y = self.knee_curve()
        result = detect_inflection(x, y, min_n=50)
        self.assertTrue(result['inflection'])
        self.assertLessEqual(abs(result['n'] - 750), 750 * 0.25)
        self.assertIsNotNone(result['kneedle'])
        self.assertIsNotNone(result['second_difference'])

    def test_straight_line_has_no_inflection(self):
        x = np.arange(1, 2001, dtype=float)
        y = 0.001 * x
        result = detect_inflection(x, y, min_n=50)
        self.assertFalse(result['inflection'])
        self.assertIsNone(result['n'])
        self.assertIsNone(result['kneedle'])
        self.assertIsNone(second_difference_peak(y, min_n=50))
        self.assertIsNone(kneedle_knee(x, y))

    def test_leading_outliers_do_not_fake_a_knee(self):
        x = np.arange(1, 2001, dtype=float)
        y = x * 0.001 + 0.4 * (x == 1)
        y = y - y.min()
        result = detect_inflection(x, y, min_n=100)
        self.assertIsNone(result['n'])

    def test_declared_criterion_names_both_tests(self):
        recipe = driver.RECIPE['inflection']
        self.assertIn('Kneedle', recipe)
        self.assertIn('second-difference', recipe)
        self.assertIn('inflection', recipe)


class FigureReplayTests(unittest.TestCase):
    def test_renderer_rejects_incomplete_csv(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'bad.csv'
            path.write_text('record,panel,panel_title,title\nfigure,,,,x\n',
                            encoding='utf-8')
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
            first = temporary / 'first.png'
            second = temporary / 'second.png'
            renderer.render(*renderer.load(path), first)
            renderer.render(*renderer.load(path), second)
            self.assertEqual(hashlib.sha256(first.read_bytes()).hexdigest(),
                             hashlib.sha256(second.read_bytes()).hexdigest())

    @staticmethod
    def synthetic_rows():
        rows = [{'record': 'figure', 'title': 'fixture'},
                {'record': 'panel', 'panel': 'a', 'panel_title': 'bars'},
                {'record': 'panel', 'panel': 'b', 'panel_title': 'hist'},
                {'record': 'panel', 'panel': 'c', 'panel_title': 'curve'},
                {'record': 'panel', 'panel': 'd', 'panel_title': 'matrix'}]
        for source, share in (('tenv', 0.1), ('gse141910', 0.2)):
            rows.append({'record': 'bar', 'source': source, 'hvg_share': share,
                         'deg_up_share': share / 2, 'deg_down_share': share / 3})
        for index in range(4):
            rows.append({'record': 'hist', 'bin_left': index, 'bin_right': index + 1,
                         'count': 10 + index})
        rows.append({'record': 'marker', 'marker_n': '', 'marker_score': '',
                     'note': 'no inflection'})
        for k in (1000, 2000):
            for n in (1, 100, 1000):
                rows.append({'record': 'curve', 'k': k, 'n': n,
                             'non_hvg_share': 0.5 + k / 10000 + n / 10000})
        for left in ('tenv_hvg_2000', 'gse141910_deg'):
            for right in ('tenv_hvg_2000', 'gse141910_deg'):
                rows.append({'record': 'matrix', 'row_source': left, 'col_source': right,
                             'jaccard': 1.0 if left == right else 0.25})
        return rows

    def test_committed_figure_replays_from_its_csv_alone(self):
        csv_path = STEP / 'figure/data/gene_axis_source_summary.csv'
        png_path = STEP / 'figure/gene_axis_source_summary.png'
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            copied = temporary / 'figure.csv'
            copied.write_bytes(csv_path.read_bytes())
            rendered = temporary / 'rendered.png'
            renderer.render(*renderer.load(copied), rendered)
            self.assertEqual(hashlib.sha256(rendered.read_bytes()).hexdigest(),
                             hashlib.sha256(png_path.read_bytes()).hexdigest())
        source = Path(renderer.__file__).read_text(encoding='utf-8')
        self.assertNotIn('read_csv', source)
        self.assertNotIn('read_excel', source)


class LeakageTests(unittest.TestCase):
    def test_guard_rejects_truth_and_heliste_paths(self):
        # The Heliste path is deliberately non-existent (the guard rejects by token,
        # not by file), so it is split to keep script/check_paths.sh from flagging a
        # literal that cannot resolve.
        for path in ('stage2/01_bases/output/tenv-truth/manifest.json',
                     'stage2/01_bases/output/tenv-truth-successor/labels.csv',
                     'pheno_data/Heliste-2021/source/' + 'labels.csv'):
            with self.assertRaises(ValueError):
                driver.forbidden_input(Path(path))
        driver.forbidden_input(Path('stage2/01_bases/output/tenv-hvg/descriptive.csv'))

    def test_poisoned_open_blocks_truth_but_not_evidence(self):
        real_open = builtins.open
        truth_reads = []

        def poison(file, *args, **kwargs):
            name = str(file).lower()
            if 'tenv-truth' in name or 'heliste' in name:
                truth_reads.append(str(file))
                raise AssertionError(f'truth read attempted: {file}')
            return real_open(file, *args, **kwargs)

        with mock.patch('builtins.open', poison):
            manifest = json.loads((STEP / 'output/gene-axis/manifest.json').read_text())
            self.assertTrue(manifest['inputs'])
            for path in manifest['inputs']:
                self.assertNotIn('tenv-truth', path.lower())
                self.assertNotIn('heliste', path.lower())
            with self.assertRaises(AssertionError):
                open('stage2/01_bases/output/tenv-truth/manifest.json')
            with self.assertRaises(AssertionError):
                pd.read_csv('stage2/01_bases/output/tenv-truth/labels.csv')
        self.assertEqual(len(truth_reads), 2)

    def test_manifests_contain_no_absolute_paths(self):
        manifest = json.loads((STEP / 'output/gene-axis/manifest.json').read_text())
        for path in manifest['inputs']:
            self.assertFalse(Path(path).is_absolute(), path)
        self.assertNotIn(str(driver.repo_root()), json.dumps(manifest['inputs']))


class SourceAdapterTests(unittest.TestCase):
    @staticmethod
    def dtoxs_fixture(path):
        header = ['CTRL.w1', 'DOX.w2', 'DOX.w3',
                  'CTRL.w1.Norm', 'DOX.w2.Norm', 'DOX.w3.Norm',
                  'logFC', 'logCPM', 'PValue', 'FDR', 'Regulation']
        table = ['\t'.join(header),
                 '\t'.join(['A', '10', '20', '30', '1', '2', '3',
                            '1.5', '2.0', '0.01', '0.02', '1']),
                 '\t'.join(['B', '10', '10', '10', '1', '1', '1',
                            '-0.5', '1.0', '0.20', '0.30', '-1'])]
        with tarfile.open(path, 'w:gz') as archive:
            payload = ('\n'.join(table) + '\n').encode()
            info = tarfile.TarInfo('Human-X-CM-Hour-48-Plate-0-Calc.CTRL-DOX.tsv')
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))

    def test_dtoxs_adapter_reads_upstream_statistics(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            path = temporary / 'deg.tar.gz'
            self.dtoxs_fixture(path)
            frame, provenance = driver.sources.dtoxs_contrasts(path)
            self.assertEqual(provenance['contrasts'], [Path('Human-X-CM-Hour-48-Plate-0-Calc.CTRL-DOX.tsv').stem])
            self.assertEqual(provenance['n_distinct_wells'], 3)
            self.assertEqual(frame.shape[0], 2)
            self.assertEqual(frame.n_control.iloc[0], 1)
            self.assertEqual(frame.n_treated.iloc[0], 2)
            self.assertAlmostEqual(frame.padj.iloc[0], 0.02)
            self.assertAlmostEqual(frame.log2fc.iloc[1], -0.5)


if __name__ == '__main__':
    unittest.main()
