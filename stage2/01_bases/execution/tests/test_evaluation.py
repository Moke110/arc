"""Shared Stage 2 evaluation contract: measured effects, direct, rank table, @d family."""
import hashlib
import unittest
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import evaluation
from bases.evaluation import (
    measured_effect_gate, direct_metrics, predicted_rank_table, threshold_grid,
    precision_at_d, recall_at_d, f1_at_d, spearman_at_d, reliability_at_d,
    reliability_curve, aggregate_direct, aggregate_reliability)

PADJ = np.array([0.01, 0.02, 0.049, 0.05, 0.01, 0.2])
LOG2FC = np.array([3.0, -2.0, 2.5, 10.0, -1.0, -5.0])


class MeasuredEffectTests(unittest.TestCase):
    def test_padj_selects_and_log2fc_supplies_direction_without_a_magnitude_gate(self):
        gate = measured_effect_gate(PADJ, LOG2FC)
        self.assertEqual(gate.signed.tolist(), [1, -1, 1, 0, -1, 0])
        self.assertEqual(gate.indices.tolist(), [0, 1, 2, 4])
        self.assertEqual(gate.count, 4)
        np.testing.assert_allclose(gate.magnitude, [3.0, 2.0, 2.5, 0.0, 1.0, 0.0])

    def test_only_padj_below_threshold_selects_and_nan_padj_is_not_measured(self):
        padj = np.array([0.05, 0.049, np.nan, 0.0])
        log2fc = np.array([-40.0, 1.999, 5.0, 2.0])
        gate = measured_effect_gate(padj, log2fc)
        self.assertEqual(gate.signed.tolist(), [0, 1, 0, 1])
        np.testing.assert_allclose(gate.magnitude, [0.0, 1.999, 0.0, 2.0])

    def test_gate_rejects_inconsistent_or_degenerate_statistics(self):
        with self.assertRaisesRegex(ValueError, 'adjusted p-values'):
            measured_effect_gate(np.array([1.2, 0.01]), LOG2FC[:2])
        with self.assertRaisesRegex(ValueError, 'direction'):
            measured_effect_gate(np.array([0.01, 0.02]), np.array([np.nan, -3.0]))
        with self.assertRaisesRegex(ValueError, 'direction'):
            measured_effect_gate(np.array([0.01, 0.02]), np.array([0.0, -3.0]))
        with self.assertRaisesRegex(ValueError, 'gene axis'):
            measured_effect_gate(np.array([0.01, 0.02]), LOG2FC[:3])


class DirectMetricTests(unittest.TestCase):
    def test_per_well_direct_metrics_match_hand_computed_values(self):
        truth = np.array([1., 2., 3., 4., 5.])
        prediction = np.array([5., 3., 4., 1., 2.])
        metrics = direct_metrics(truth, prediction)
        self.assertEqual(tuple(metrics.values), ('PCC', 'Spearman', 'RMSE', 'R2'))
        self.assertEqual(metrics.reasons, {})
        self.assertAlmostEqual(metrics.values['PCC'], -0.8)
        self.assertAlmostEqual(metrics.values['Spearman'], -0.8)
        self.assertAlmostEqual(metrics.values['RMSE'], np.sqrt(7.2))
        self.assertAlmostEqual(metrics.values['R2'], -2.6)

    def test_constant_prediction_and_constant_truth_are_nan_with_distinct_reasons(self):
        truth = np.array([1., 2., 3., 4., 5.])
        constant_prediction = direct_metrics(truth, np.full(5, 3.))
        self.assertTrue(np.isnan(constant_prediction.values['PCC']))
        self.assertTrue(np.isnan(constant_prediction.values['Spearman']))
        self.assertTrue(np.isnan(constant_prediction.values['R2']))
        self.assertEqual(constant_prediction.reasons['PCC'], 'constant_prediction')
        self.assertEqual(constant_prediction.reasons['Spearman'], 'constant_prediction')
        self.assertEqual(constant_prediction.reasons['R2'], 'constant_prediction')
        self.assertAlmostEqual(constant_prediction.values['RMSE'], np.sqrt(2.))
        self.assertNotIn('RMSE', constant_prediction.reasons)

        constant_truth = direct_metrics(np.full(5, 2.), np.array([1., 2., 3., 4., 5.]))
        self.assertTrue(np.isnan(constant_truth.values['PCC']))
        self.assertTrue(np.isnan(constant_truth.values['R2']))
        self.assertEqual(constant_truth.reasons['PCC'], 'constant_truth')
        self.assertEqual(constant_truth.reasons['R2'], 'constant_truth')
        self.assertAlmostEqual(constant_truth.values['RMSE'], np.sqrt(3.))

    def test_direct_metrics_require_one_finite_axis(self):
        with self.assertRaisesRegex(ValueError, 'gene axis'):
            direct_metrics(np.ones(3), np.ones(4))
        with self.assertRaisesRegex(ValueError, 'finite'):
            direct_metrics(np.array([1., np.nan, 3.]), np.ones(3))


class PredictedRankTableTests(unittest.TestCase):
    def test_absolute_rank_table_orders_by_abs_desc_and_breaks_ties_by_gene_index(self):
        table = predicted_rank_table(np.array([10., 6., -6., 1., 3.]))
        self.assertEqual(table.gene_index.tolist(), [0, 1, 2, 4, 3])
        self.assertEqual(table.rank.tolist(), [1, 2, 3, 4, 5])
        self.assertEqual(table.direction.tolist(), [1, 1, -1, 1, 1])
        self.assertEqual(table.absolute_delta.tolist(), [10., 6., 6., 3., 1.])

    def test_rank_table_covers_the_whole_axis_and_keeps_zero_predictions_last(self):
        table = predicted_rank_table(np.array([0., -2., 0., 5.]))
        self.assertEqual(table.gene_index.tolist(), [3, 1, 0, 2])
        self.assertEqual(table.direction.tolist(), [1, -1, 0, 0])
        self.assertEqual(sorted(table.rank.tolist()), [1, 2, 3, 4])

    def test_rank_table_rejects_malformed_predictions(self):
        with self.assertRaisesRegex(ValueError, 'axis'):
            predicted_rank_table(np.zeros((2, 2)))
        with self.assertRaisesRegex(ValueError, 'finite'):
            predicted_rank_table(np.array([1., np.inf]))


def _gate(padj, log2fc):
    return measured_effect_gate(np.asarray(padj, dtype=float), np.asarray(log2fc, dtype=float))


def _curve(delta, padj, log2fc, max_absolute_delta=10.0):
    return reliability_curve(np.asarray(delta, dtype=float), _gate(padj, log2fc),
                             threshold_grid(max_absolute_delta))


class ThresholdGridTests(unittest.TestCase):
    def test_grid_is_one_hundred_scaled_points_with_open_ends(self):
        grid = threshold_grid(8.0)
        self.assertEqual(len(grid.d_prime), 100)
        np.testing.assert_array_equal(grid.d_prime, np.arange(1, 101) / 100)
        self.assertEqual(grid.d_prime[0], 0.01)
        self.assertEqual(grid.d_prime[-1], 1.0)
        self.assertAlmostEqual(grid.thresholds[0], 0.08)
        self.assertAlmostEqual(grid.thresholds[-1], 8.0)
        np.testing.assert_allclose(grid.thresholds, np.arange(1, 101) * 0.01 * 8.0)
        self.assertEqual(grid.max_absolute_delta, 8.0)

    def test_two_runs_with_different_m_share_one_d_prime_axis(self):
        small, large = threshold_grid(3.0), threshold_grid(0.4)
        np.testing.assert_array_equal(small.d_prime, large.d_prime)
        self.assertFalse(np.allclose(small.thresholds, large.thresholds))

    def test_grid_rejects_nonpositive_or_nonfinite_normalisers(self):
        for value in (0.0, -1.0, np.nan, np.inf):
            with self.assertRaisesRegex(ValueError, 'max_absolute_delta'):
                threshold_grid(value)


class ReliabilityMetricTests(unittest.TestCase):
    def test_four_grid_cells_take_their_declared_totals_and_stay_finite(self):
        grid = threshold_grid(10.0)
        empty_none = reliability_curve(np.zeros(3), _gate([.2] * 3, [5.] * 3), grid)
        empty_measured = reliability_curve(np.zeros(3), _gate([.01, .2, .2], [1., 5., 5.]), grid)
        window_none = reliability_curve(np.array([10., 0., 0.]), _gate([.2] * 3, [5.] * 3), grid)
        window_measured = reliability_curve(
            np.array([10., 5., 0.]), _gate([.01, .2, .01], [2., 5., 1.]), grid)

        self.assertEqual(empty_none.precision[99], 1.0)
        self.assertEqual(empty_none.recall[99], 1.0)
        self.assertEqual(empty_none.f1[99], 1.0)
        self.assertEqual(empty_none.spearman[99], 1.0)
        self.assertEqual(empty_none.reliability[99], 1.0)

        for curve in (empty_measured, window_none):
            self.assertEqual(curve.precision[99], 0.0)
            self.assertEqual(curve.recall[99], 0.0)
            self.assertEqual(curve.f1[99], 0.0)
            self.assertEqual(curve.reliability[99], 0.0)

        self.assertEqual(window_measured.n_window[39], 2)
        self.assertEqual(window_measured.n_hits[39], 1)
        self.assertEqual(window_measured.precision[39], 0.5)
        self.assertEqual(window_measured.recall[39], 0.5)
        self.assertEqual(window_measured.f1[39], 0.5)
        self.assertEqual(window_measured.spearman[39], 1.0)
        self.assertEqual(window_measured.reliability[39], 0.5)

        for curve in (empty_none, empty_measured, window_none, window_measured):
            self.assertFalse(hasattr(curve, 'reasons'))
            for attribute in ('precision', 'recall', 'f1', 'spearman', 'reliability'):
                self.assertTrue(np.isfinite(getattr(curve, attribute)).all())

    def test_perfect_table_scores_reliability_one(self):
        curve = _curve([4., 3., 2., 0.], [.01, .01, .01, .2], [4., 3., 2., 5.], 4.0)
        index = 49  # d' = 0.50, so d = 2.0
        self.assertEqual(curve.n_window[index], 3)
        self.assertEqual(curve.n_hits[index], 3)
        self.assertEqual(curve.precision[index], 1.0)
        self.assertEqual(curve.recall[index], 1.0)
        self.assertEqual(curve.f1[index], 1.0)
        self.assertEqual(curve.spearman[index], 1.0)
        self.assertEqual(curve.reliability[index], 1.0)

    def test_direction_inverted_perfect_set_scores_zero(self):
        curve = _curve([-4., -3., -2., 0.], [.01, .01, .01, .2], [4., 3., 2., 5.], 4.0)
        index = 49  # d' = 0.50, so d = 2.0 and the full measured set is still emitted
        self.assertEqual(curve.n_window[index], 3)
        self.assertEqual(curve.n_hits[index], 0)
        self.assertEqual(curve.f1[index], 0.0)
        self.assertEqual(curve.reliability[index], 0.0)

    def test_f1_is_harmonic_so_one_hit_in_three_is_half_not_two_thirds(self):
        curve = _curve([10., 0., 0., 0.], [.01, .01, .01, .2], [3., 2., 1., 5.])
        index = 99  # d = 10.0 selects only the strongest gene
        self.assertEqual(curve.precision[index], 1.0)
        self.assertAlmostEqual(curve.recall[index], 1 / 3)
        self.assertAlmostEqual(curve.f1[index], 0.5)
        self.assertNotAlmostEqual(curve.f1[index], (1.0 + 1 / 3) / 2)
        self.assertAlmostEqual(curve.reliability[index], 0.5)

    def test_spearman_is_neutral_below_three_hits_and_hand_computed_at_three(self):
        self.assertEqual(spearman_at_d(np.array([]), np.array([])), 1.0)
        self.assertEqual(spearman_at_d(np.array([3.]), np.array([1.])), 1.0)
        self.assertEqual(spearman_at_d(np.array([3., 2.]), np.array([1., 2.])), 1.0)
        # ranks 1,2,3 against 1,3,2: 1 - 6*2/(3*8) = 0.5.
        self.assertAlmostEqual(spearman_at_d(np.array([1., 2., 3.]), np.array([1., 3., 2.])), 0.5)
        self.assertAlmostEqual(reliability_at_d(0.6, 0.5), 0.45)
        # No rank variation is uninformative, so the term stays neutral and finite.
        self.assertEqual(spearman_at_d(np.array([2., 2., 2.]), np.array([1., 2., 3.])), 1.0)

    def test_correct_set_in_reverse_magnitude_order_loses_the_ordering_term(self):
        # The window is the measured set, but predicted magnitudes rise as measured fall.
        curve = _curve([2., 3., 4., 0.], [.01, .01, .01, .2], [4., 3., 2., 5.], 4.0)
        self.assertEqual(curve.f1[49], 1.0)
        self.assertAlmostEqual(curve.spearman[49], -1.0)
        self.assertEqual(curve.reliability[49], 0.0)

    def test_strictest_grid_point_selects_the_gene_that_defines_m(self):
        curve = _curve([25., 10., 0.], [.2, .2, .2], [5., 5., 5.], 25.0)
        self.assertAlmostEqual(curve.thresholds[-1], 25.0)
        self.assertEqual(curve.n_window[-1], 1)
        self.assertEqual(curve.n_hits[-1], 0)


class AggregationTests(unittest.TestCase):
    def test_direct_run_mean_is_equal_weight_over_defined_wells(self):
        first = direct_metrics(np.array([1., 2., 3., 4., 5.]), np.array([5., 3., 4., 1., 2.]))
        second = direct_metrics(np.array([1., 2., 3., 4., 5.]), np.array([1., 2., 3., 4., 5.]))
        undefined = direct_metrics(np.array([1., 2., 3., 4., 5.]), np.full(5, 3.))
        aggregate = aggregate_direct([first, second, undefined])
        self.assertEqual(aggregate.counts['wells'], 3)
        self.assertEqual(aggregate.counts['PCC'], 2)
        self.assertEqual(aggregate.counts['RMSE'], 3)
        self.assertAlmostEqual(aggregate.values['PCC'], 0.1)
        self.assertAlmostEqual(aggregate.values['Spearman'], 0.1)
        self.assertAlmostEqual(aggregate.values['R2'], -0.8)
        self.assertAlmostEqual(aggregate.values['RMSE'], (np.sqrt(7.2) + np.sqrt(2.)) / 3)
        self.assertEqual(aggregate.reasons, {})

    def test_direct_aggregation_reports_undefined_metrics_with_a_reason(self):
        wells = [direct_metrics(np.full(3, 2.), np.array([1., 2., 3.]))]
        aggregate = aggregate_direct(wells)
        self.assertTrue(np.isnan(aggregate.values['PCC']))
        self.assertEqual(aggregate.counts['PCC'], 0)
        self.assertEqual(aggregate.reasons['PCC'], 'constant_truth')
        self.assertEqual(aggregate.counts['RMSE'], 1)

    def test_macro_is_equal_weight_over_every_passed_test_group_with_a_null_split(self):
        changed = _curve([4., 3., 2., 0.], [.01, .01, .01, .2], [4., 3., 2., 5.])
        null_empty = _curve([0., 0., 0., 0.], [.2] * 4, [5.] * 4)
        null_hallucinated = _curve([10., 0., 0., 0.], [.2] * 4, [5.] * 4)
        aggregate = aggregate_reliability([changed, null_empty, null_hallucinated])

        self.assertEqual(aggregate.counts, {'groups': 3, 'null_groups': 2, 'changed_groups': 1})
        np.testing.assert_allclose(
            aggregate.reliability,
            (changed.reliability + null_empty.reliability + null_hallucinated.reliability) / 3)
        np.testing.assert_allclose(
            aggregate.null_reliability,
            (null_empty.reliability + null_hallucinated.reliability) / 2)
        np.testing.assert_allclose(aggregate.changed_reliability, changed.reliability)
        self.assertAlmostEqual(aggregate.reliability[99], 1 / 3)
        np.testing.assert_array_equal(aggregate.d_prime, changed.grid.d_prime)
        for attribute in ('precision', 'recall', 'f1', 'spearman', 'reliability',
                          'null_reliability', 'changed_reliability'):
            self.assertTrue(np.isfinite(getattr(aggregate, attribute)).all())
        self.assertFalse(hasattr(aggregate, 'reasons'))

    def test_macro_uses_only_the_passed_groups(self):
        first = _curve([4., 3., 2., 0.], [.01, .01, .01, .2], [4., 3., 2., 5.])
        second = _curve([0., 0., 0., 0.], [.2] * 4, [5.] * 4)
        pair = aggregate_reliability([first, second])
        self.assertEqual(pair.counts['groups'], 2)
        np.testing.assert_allclose(pair.reliability, (first.reliability + second.reliability) / 2)

    def test_macro_requires_at_least_one_group_and_one_shared_grid(self):
        with self.assertRaisesRegex(ValueError, 'at least one'):
            aggregate_reliability([])
        first = _curve([4., 3., 2., 0.], [.01, .01, .01, .2], [4., 3., 2., 5.], 10.0)
        other_run = _curve([4., 3., 2., 0.], [.01, .01, .01, .2], [4., 3., 2., 5.], 20.0)
        with self.assertRaisesRegex(ValueError, 'threshold grid'):
            aggregate_reliability([first, other_run])


class RetiredVocabularyTests(unittest.TestCase):
    """The indirect family and its NaN-plus-reason machinery are gone (ADR 0064).

    The retired tokens below are the guard itself; nothing else in this file may
    use them.
    """

    RETIRED_TOKENS = (
        'Precision@k', 'Spearman@k', 'max_k@P80S40', 'Recall@max_k',
        'INDIRECT_METRICS', 'EFFECT_THRESHOLD', 'K_MIN_FOR_SPEARMAN', 'K_MAX',
        'PRECISION_TARGET', 'SPEARMAN_TARGET', 'empty_gate', 'k_below_3',
        'no_qualifying_k', 'no_defined_groups')

    def test_no_retired_name_gate_or_reason_remains_in_the_implementation(self):
        source = Path(evaluation.__file__).read_text()
        for token in self.RETIRED_TOKENS:
            self.assertNotIn(token, source)
        for name in ('truth_gate', 'GatedTruth', 'GroupIndirectMetrics', 'IndirectAggregate',
                     'indirect_metrics', 'aggregate_indirect'):
            self.assertFalse(hasattr(evaluation, name), name)

    def test_the_new_family_and_its_docstring_are_declared(self):
        self.assertEqual(evaluation.RELIABILITY_METRICS,
                         ('Precision@d', 'Recall@d', 'F1@d', 'Spearman@d', 'Reliability@d'))
        self.assertEqual(evaluation.MIN_HITS_FOR_SPEARMAN, 3)
        self.assertEqual(evaluation.SPEARMAN_NEUTRAL, 1.0)
        self.assertIn('ADR 0064', evaluation.__doc__)


class RetiredAuprcTests(unittest.TestCase):
    """The retired AUPRC/PR implementation stays byte-identical and importable."""

    def test_retired_module_is_untouched_and_still_scores(self):
        import bases.metrics as retired

        digest = hashlib.sha256(Path(retired.__file__).read_bytes()).hexdigest()
        self.assertEqual(digest, 'f9278dc6aa31e3a2d32a534303d4a8051bc6ec89fdd858b6ec9c49960ee541d8')
        self.assertEqual(retired.METRICS, ('AUPRC', 'AUPRC_positive', 'AUPRC_negative'))
        curve = retired.average_precision(np.array([2., -2., 1.]), np.array([1, 0, 1]), 4)
        self.assertAlmostEqual(curve.ap, .5 / 4 + (2 / 3) / 4)
        continuous = retired.continuous(np.array([1., 2., 3.]), np.array([2., 4., 6.]))
        self.assertEqual(continuous['PCC'], 1.)
        self.assertEqual(continuous['RMSE'], np.sqrt(14 / 3))
