"""Shared Stage 2 evaluation contract: gate, direct, rank table, indirect, aggregation."""
import hashlib
import unittest
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.evaluation import (
    truth_gate, direct_metrics, predicted_rank_table, indirect_metrics,
    aggregate_direct, aggregate_indirect)

PADJ = np.array([0.01, 0.02, 0.049, 0.05, 0.01, 0.2])
LOG2FC = np.array([3.0, -2.0, 2.5, 10.0, 1.0, -5.0])
DELTA_MEAN = np.array([5.0, -3.0, 2.0, 4.0, 1.0, -1.0])


class TruthGateTests(unittest.TestCase):
    def test_gate_requires_padj_and_fold_change_and_takes_direction_from_delta_mean(self):
        gate = truth_gate(PADJ, LOG2FC, DELTA_MEAN)
        self.assertEqual(gate.signed.tolist(), [1, -1, 1, 0, 0, 0])
        self.assertEqual(gate.indices.tolist(), [0, 1, 2])
        self.assertEqual(gate.count, 3)

    def test_gate_ignores_non_significant_effect_and_threshold_edges(self):
        padj = np.array([0.05, 0.049, 0.01, 0.2])
        log2fc = np.array([40.0, 1.999, -2.0, 7.0])
        delta_mean = np.array([8.0, 4.0, -6.0, 3.0])
        gate = truth_gate(padj, log2fc, delta_mean)
        self.assertEqual(gate.signed.tolist(), [0, 0, -1, 0])

    def test_gate_rejects_inconsistent_or_degenerate_statistics(self):
        with self.assertRaisesRegex(ValueError, 'adjusted p-values'):
            truth_gate(np.array([1.2, 0.01]), LOG2FC[:2], DELTA_MEAN[:2])
        with self.assertRaisesRegex(ValueError, 'significant gene'):
            truth_gate(np.array([0.01, 0.02]), np.array([np.nan, -3.0]), DELTA_MEAN[:2])
        with self.assertRaisesRegex(ValueError, 'measured delta'):
            truth_gate(np.array([0.01, 0.02]), np.array([3.0, -3.0]), np.array([0.0, -3.0]))


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


def _gated(padj, log2fc, delta_mean):
    return truth_gate(np.asarray(padj, dtype=float), np.asarray(log2fc, dtype=float),
                      np.asarray(delta_mean, dtype=float))


class IndirectMetricTests(unittest.TestCase):
    def test_exact_k_curve_with_tie_in_absolute_predicted_delta(self):
        # genes 1 and 2 tie at |Δ̂| = 6, so the table orders them by gene index.
        true_delta = np.array([5., -4., 3., 2.5])
        gated = _gated([.01, .01, .01, .2], [3., -2., 3., 5.], true_delta)
        table = predicted_rank_table(np.array([10., 6., -6., 1.]))
        self.assertEqual(table.gene_index.tolist(), [0, 1, 2, 3])
        metrics = indirect_metrics(table, gated, true_delta, k_max=4)
        self.assertEqual(metrics.k.tolist(), [1, 2, 3, 4])
        self.assertEqual(metrics.precision[0], 1.)
        self.assertAlmostEqual(metrics.precision[1], 1 / 2)
        self.assertAlmostEqual(metrics.precision[2], 1 / 3)
        self.assertAlmostEqual(metrics.precision[3], 1 / 4)
        self.assertTrue(np.isnan(metrics.spearman[0]))
        self.assertTrue(np.isnan(metrics.spearman[1]))
        self.assertAlmostEqual(metrics.spearman[2], 0.5)
        self.assertAlmostEqual(metrics.spearman[3], 0.2)
        self.assertEqual(metrics.reasons['Spearman@k'], 'k_below_3')
        self.assertTrue(np.isnan(metrics.max_k))
        self.assertTrue(np.isnan(metrics.recall))
        self.assertEqual(metrics.reasons['max_k@P80S40'], 'no_qualifying_k')
        self.assertEqual(metrics.reasons['Recall@max_k'], 'no_qualifying_k')
        self.assertEqual(metrics.matched, 1)

    def test_max_k_and_recall_on_a_qualifying_group(self):
        # Absolute predicted order: 0,1,2,3,10,4,11,5,... ; gated truth 0..5.
        predicted = np.array([12., 11., 10., 9., 7., 5., 4., 3., 2., 1., 8., 6.])
        true_delta = np.array([6., 5., 4., 3., 2., 1., 0., 0., 0., 0., 0., 0.])
        padj = [.01] * 6 + [.2] * 6
        log2fc = [3., 2., 2.5, 3., 2., 2.] + [5.] * 6
        gated = _gated(padj, log2fc, true_delta)
        metrics = indirect_metrics(predicted_rank_table(predicted), gated, true_delta, k_max=12)
        self.assertEqual(metrics.precision[:4].tolist(), [1., 1., 1., 1.])
        self.assertAlmostEqual(metrics.precision[4], .8)
        self.assertAlmostEqual(metrics.precision[5], 5 / 6)
        self.assertAlmostEqual(metrics.spearman[2], 1.)
        self.assertAlmostEqual(metrics.spearman[3], 1.)
        self.assertAlmostEqual(metrics.spearman[4], 1.)
        self.assertAlmostEqual(metrics.spearman[5], 33 / 35)
        self.assertEqual(metrics.max_k, 6.)
        self.assertAlmostEqual(metrics.recall, 5 / 6)
        self.assertEqual(metrics.matched, 6)
        self.assertNotIn('max_k@P80S40', metrics.reasons)

    def test_full_k_range_reaches_one_hundred(self):
        count = 120
        predicted = np.arange(count, 0, -1, dtype=float)
        true_delta = np.arange(count, 0, -1, dtype=float)
        padj = np.full(count, .2)
        padj[:60] = .01
        padj[110] = .01
        log2fc = np.full(count, 5.)
        gated = _gated(padj, log2fc, true_delta)
        metrics = indirect_metrics(predicted_rank_table(predicted), gated, true_delta)
        self.assertEqual(len(metrics.k), 100)
        self.assertEqual(metrics.precision[0], 1.)
        self.assertEqual(metrics.precision[59], 1.)
        self.assertAlmostEqual(metrics.precision[60], 60 / 61)
        self.assertAlmostEqual(metrics.precision[74], 60 / 75)
        self.assertAlmostEqual(metrics.precision[75], 60 / 76)
        self.assertAlmostEqual(metrics.precision[99], .6)
        np.testing.assert_allclose(metrics.spearman[2:], np.ones(98))
        self.assertEqual(metrics.max_k, 75.)
        self.assertAlmostEqual(metrics.recall, 60 / 61)
        self.assertEqual(metrics.matched, 61)

    def test_empty_gate_is_nan_with_its_own_reason(self):
        true_delta = np.array([5., 4., -3.])
        gated = _gated([.2, .2, .2], [5., 5., 5.], true_delta)
        metrics = indirect_metrics(predicted_rank_table(np.array([3., 2., 1.])), gated, true_delta, k_max=3)
        self.assertTrue(np.isnan(metrics.precision).all())
        self.assertTrue(np.isnan(metrics.spearman).all())
        self.assertTrue(np.isnan(metrics.max_k))
        self.assertTrue(np.isnan(metrics.recall))
        for name in ('Precision@k', 'Spearman@k', 'max_k@P80S40', 'Recall@max_k'):
            self.assertEqual(metrics.reasons[name], 'empty_gate')

    def test_constant_prediction_and_constant_truth_are_distinct_reasons(self):
        true_delta = np.array([5., 4., 3., 2.])
        gated = _gated([.01] * 3 + [.2], [3.] * 4, true_delta)
        constant_prediction = indirect_metrics(
            predicted_rank_table(np.zeros(4)), gated, true_delta, k_max=4)
        self.assertTrue(np.isnan(constant_prediction.spearman).all())
        self.assertEqual(constant_prediction.reasons['Spearman@k'], 'constant_prediction')
        self.assertEqual(constant_prediction.reasons['max_k@P80S40'], 'no_qualifying_k')

        constant_truth = indirect_metrics(
            predicted_rank_table(np.array([4., 3., 2., 1.])), gated, np.full(4, 2.), k_max=4)
        self.assertTrue(np.isnan(constant_truth.spearman).all())
        self.assertEqual(constant_truth.reasons['Spearman@k'], 'constant_truth')

    def test_indirect_metrics_require_one_axis_and_available_k(self):
        gated = _gated([.01, .01], [3., 3.], np.array([2., 3.]))
        table = predicted_rank_table(np.array([2., 1.]))
        with self.assertRaisesRegex(ValueError, 'gene axis'):
            indirect_metrics(table, gated, np.array([2., 3.]), k_max=3)
        with self.assertRaisesRegex(ValueError, 'k_max'):
            indirect_metrics(table, gated, np.array([2., 3.]), k_max=0)


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

    def test_indirect_run_mean_is_equal_weight_over_defined_groups(self):
        predicted = np.array([12., 11., 10., 9., 7., 5., 4., 3., 2., 1., 8., 6.])
        true_delta = np.array([6., 5., 4., 3., 2., 1., 0., 0., 0., 0., 0., 0.])
        qualifying = indirect_metrics(
            predicted_rank_table(predicted),
            _gated([.01] * 6 + [.2] * 6, [3., 2., 2.5, 3., 2., 2.] + [5.] * 6, true_delta),
            true_delta, k_max=12)
        # One gated truth gene at the top cannot reach precision .8 for k >= 2.
        descending = np.arange(12., 0., -1)
        unqualified = indirect_metrics(
            predicted_rank_table(descending),
            _gated([.01] + [.2] * 11, [3.] + [5.] * 11, np.array([5.] + [0.] * 11)),
            np.array([5.] + [0.] * 11), k_max=12)
        empty = indirect_metrics(
            predicted_rank_table(descending),
            _gated([.2] * 12, [5.] * 12, np.arange(12., 0., -1)),
            np.arange(12., 0., -1), k_max=12)
        aggregate = aggregate_indirect([qualifying, unqualified, empty])
        self.assertEqual(aggregate.counts['groups'], 3)
        self.assertEqual(aggregate.k.tolist(), list(range(1, 13)))
        self.assertEqual(aggregate.counts['Precision@k'].tolist(), [2] * 12)
        np.testing.assert_allclose(aggregate.precision, (qualifying.precision + unqualified.precision) / 2)
        self.assertEqual(aggregate.counts['Spearman@k'][:2].tolist(), [0, 0])
        self.assertEqual(aggregate.counts['Spearman@k'][2], 2)
        self.assertAlmostEqual(aggregate.max_k, 6.)
        self.assertEqual(aggregate.counts['max_k@P80S40'], 1)
        self.assertAlmostEqual(aggregate.recall, 5 / 6)
        self.assertEqual(aggregate.counts['Recall@max_k'], 1)
        self.assertEqual(aggregate.reasons, {})

    def test_indirect_aggregation_with_every_group_undefined_reason(self):
        empty = indirect_metrics(
            predicted_rank_table(np.array([3., 2., 1.])),
            _gated([.2, .2, .2], [5., 5., 5.], np.array([5., 4., 3.])),
            np.array([5., 4., 3.]), k_max=3)
        aggregate = aggregate_indirect([empty])
        self.assertTrue(np.isnan(aggregate.precision).all())
        self.assertTrue(np.isnan(aggregate.spearman).all())
        self.assertTrue(np.isnan(aggregate.max_k))
        self.assertTrue(np.isnan(aggregate.recall))
        self.assertEqual(aggregate.counts['Precision@k'].tolist(), [0, 0, 0])
        self.assertEqual(aggregate.counts['max_k@P80S40'], 0)
        for name in ('Precision@k', 'Spearman@k', 'max_k@P80S40', 'Recall@max_k'):
            self.assertEqual(aggregate.reasons[name], 'no_defined_groups')

    def test_indirect_aggregation_requires_matching_k_axes(self):
        empty = indirect_metrics(
            predicted_rank_table(np.array([3., 2., 1.])),
            _gated([.2, .2, .2], [5., 5., 5.], np.array([5., 4., 3.])), np.array([5., 4., 3.]), k_max=3)
        other = indirect_metrics(
            predicted_rank_table(np.array([4., 3., 2., 1.])),
            _gated([.2] * 4, [5.] * 4, np.array([5., 4., 3., 2.])),
            np.array([5., 4., 3., 2.]), k_max=4)
        with self.assertRaisesRegex(ValueError, 'k axis'):
            aggregate_indirect([empty, other])


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
