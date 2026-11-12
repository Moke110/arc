"""Fail-closed run-set and coverage contract of the tenv evaluation driver (#164)."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

EXECUTION = Path(__file__).resolve().parents[1]
SCRIPT = EXECUTION / 'script'
sys.path.insert(0, str(SCRIPT))
from evaluation_driver import (  # noqa: E402
    DriverError, figure_entries, learning_models, load_run_set, per_group_rows,
    per_run_rows, per_well_rows, score_run, validate_run)

GENES = 18164
CASES = ['tenv:temp_env::Plate1-A01', 'tenv:temp_env::Plate1-A02',
         'tenv:temp_env::Plate1-A03']


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def write_run(root, entry, fold, seed, *, cases=CASES, prediction=None,
              status='complete', history=None, identity=None, write_terminal=True,
              write_cases=True, write_prediction=True):
    run = Path(root) / entry / f'fold_{fold}' / f'seed_{seed}'
    run.mkdir(parents=True, exist_ok=True)
    if prediction is None:
        prediction = np.arange(len(cases) * GENES, dtype=np.float64).reshape(len(cases), GENES)
    if write_prediction:
        np.save(run / 'prediction.npy', np.asarray(prediction))
    if write_cases:
        write_json(run / 'case_ids.json', list(cases))
    if write_terminal:
        write_json(run / 'terminal.json', {'status': status, 'prediction': 'prediction.npy'})
    write_json(run / 'identity.json', {
        'identity': identity or {'model': entry, 'fold': fold, 'fold_row': fold,
                                 'seed': seed, 'seed_row': 0}})
    if history is not None:
        write_json(run / 'checkpoint.json', {'history': history})
    return run


class RunSetContract(unittest.TestCase):
    def test_complete_run_is_accepted_and_reordered_to_the_expected_case_axis(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expected = [CASES[2], CASES[0], CASES[1]]
            prediction = np.stack([np.full(GENES, value) for value in (0., 1., 2.)])
            write_run(root, 'prnet', 0, 17, cases=list(reversed(expected)), prediction=prediction)
            runs = load_run_set(root, ['prnet'], [0], [17],
                                {0: expected}, GENES)
            self.assertEqual(list(runs), [('prnet', 0, 17)])
            run = runs[('prnet', 0, 17)]
            self.assertEqual(run.case_ids, expected)
            np.testing.assert_array_equal(run.prediction[:, 0], [2., 1., 0.])

    def test_missing_entry_of_the_run_set_fails_completion(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_run(root, 'prnet', 0, 17)
            with self.assertRaisesRegex(DriverError, 'incomplete run set'):
                load_run_set(root, ['prnet', 'mean'], [0], [17], {0: CASES}, GENES)

    def test_partially_written_run_fails_completion(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_run(root, 'prnet', 0, 17, write_terminal=False)
            with self.assertRaisesRegex(DriverError, 'partially written'):
                load_run_set(root, ['prnet'], [0], [17], {0: CASES}, GENES)
            with self.assertRaisesRegex(DriverError, 'status'):
                write_run(root, 'prnet', 0, 17, status='failed')
                load_run_set(root, ['prnet'], [0], [17], {0: CASES}, GENES)

    def test_run_without_prediction_or_case_identity_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_run(root, 'prnet', 0, 17, write_prediction=False)
            with self.assertRaisesRegex(DriverError, 'prediction'):
                validate_run(root / 'prnet/fold_0/seed_17', CASES, GENES)
            write_run(root, 'prnet', 0, 18, write_cases=False)
            with self.assertRaisesRegex(DriverError, 'case_ids'):
                validate_run(root / 'prnet/fold_0/seed_18', CASES, GENES)

    def test_prediction_set_must_cover_every_fold_case(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_run(root, 'prnet', 0, 17, cases=CASES[:2])
            with self.assertRaisesRegex(DriverError, 'coverage'):
                validate_run(root / 'prnet/fold_0/seed_17', CASES, GENES)
            write_run(root, 'prnet', 0, 17,
                      cases=[*CASES[:2], 'tenv:temp_env::Plate1-Z99', CASES[2]])
            with self.assertRaisesRegex(DriverError, 'coverage'):
                validate_run(root / 'prnet/fold_0/seed_17', CASES, GENES)


def _closed_forms():
    """Hand-derived second moments of the reversed-order fixture, no implementation involved."""
    n = GENES
    first = n * (n * n - 1) / 12.0                 # SST of 1..n
    m = n - 1
    sse = n * m * m - 4 * m * (n * (n - 1) / 2) + 4 * (n * (n - 1) * (2 * n - 1) / 6)
    return first, sse


def expected_reversed_rmse():
    return float(np.sqrt(_closed_forms()[1] / GENES))


def expected_reversed_r2():
    first, sse = _closed_forms()
    return 1.0 - sse / first


class TruthRows:
    """Minimal fixture of the accepted truth rows the driver scores against."""

    def __init__(self, rows):
        self.rows = rows

    def row(self, group_id):
        return self.rows[group_id]


def synthetic_fold():
    """A known-value 120-gene fold: two test wells in two parent-dose groups."""
    axis = np.arange(GENES, dtype=np.float64)
    measured = np.stack([GENES - axis, axis + 1.0])       # well A descends, well B ascends
    predicted = np.stack([GENES - axis, GENES - axis])    # both predict the descending order
    padj = np.full((2, GENES), 0.5)
    log2fc = np.zeros((2, GENES))
    padj[0, :60] = 1e-8
    padj[0, 110] = 1e-8
    log2fc[0, :60] = 3.0
    log2fc[0, 110] = 3.0
    padj[1, 0] = 1e-8
    log2fc[1, 0] = 3.0
    rows = {
        'group-A': {'padj': padj[0], 'log2_fold_change': log2fc[0],
                    'status': 'eligible', 'reason': ''},
        'group-B': {'padj': padj[1], 'log2_fold_change': log2fc[1],
                    'status': 'eligible', 'reason': ''},
    }
    expected = [
        {'sample_id': CASES[0], 'group_id': 'group-A',
         'drug_parent': 'PARENT-A', 'dose_um': 1.0},
        {'sample_id': CASES[1], 'group_id': 'group-B',
         'drug_parent': 'PARENT-B', 'dose_um': 2.0},
    ]
    return expected, measured, predicted, TruthRows(rows)


class ScoringContract(unittest.TestCase):
    def test_direct_metrics_match_hand_computed_values_per_well(self):
        expected, measured, predicted, truth = synthetic_fold()
        score = score_run(expected, predicted, measured, truth)
        well = {row['sample_id']: row['values'] for row in score.well_results}
        first = well[CASES[0]]
        self.assertAlmostEqual(first['PCC'], 1.)
        self.assertAlmostEqual(first['Spearman'], 1.)
        self.assertEqual(first['RMSE'], 0.)
        self.assertAlmostEqual(first['R2'], 1.)
        second = well[CASES[1]]
        self.assertAlmostEqual(second['PCC'], -1.)
        self.assertAlmostEqual(second['Spearman'], -1.)
        self.assertAlmostEqual(second['R2'], expected_reversed_r2())
        self.assertAlmostEqual(second['RMSE'], expected_reversed_rmse())

    def test_indirect_metrics_match_hand_computed_values_per_group(self):
        expected, measured, predicted, truth = synthetic_fold()
        score = score_run(expected, predicted, measured, truth)
        groups = {meta['group_id']: result for result, meta
                  in zip(score.group_results, score.group_meta)}
        first = groups['group-A']
        self.assertEqual(first.precision[0], 1.)
        self.assertEqual(first.precision[59], 1.)
        self.assertAlmostEqual(first.precision[60], 60 / 61)
        self.assertAlmostEqual(first.precision[74], 60 / 75)
        self.assertAlmostEqual(first.precision[99], .6)
        self.assertAlmostEqual(first.spearman[2], 1.)
        self.assertEqual(first.max_k, 75.)
        self.assertAlmostEqual(first.recall, 60 / 61)
        self.assertEqual(first.matched, 61)
        second = groups['group-B']
        self.assertEqual(second.precision[0], 1.)
        self.assertAlmostEqual(second.precision[1], .5)
        self.assertAlmostEqual(second.spearman[2], -1.)
        self.assertTrue(np.isnan(second.max_k))
        self.assertEqual(second.reasons['max_k@P80S40'], 'no_qualifying_k')

    def test_run_tables_carry_macro_means_counts_and_reasons(self):
        expected, measured, predicted, truth = synthetic_fold()
        score = score_run(expected, predicted, measured, truth)
        run = per_run_rows([score])[0]
        self.assertEqual(run['wells'], 2)
        self.assertEqual(run['groups'], 2)
        self.assertAlmostEqual(float(run['PCC']), 0.)
        self.assertAlmostEqual(float(run['Spearman']), 0.)
        self.assertAlmostEqual(float(run['R2']), -1.)
        self.assertAlmostEqual(float(run['RMSE']), (expected_reversed_rmse()) / 2)
        self.assertEqual(int(run['PCC_defined']), 2)
        self.assertEqual(int(run['PCC_undefined']), 0)
        self.assertAlmostEqual(float(run['max_k@P80S40']), 75.)
        self.assertEqual(int(run['max_k@P80S40_defined']), 1)
        self.assertEqual(int(run['max_k@P80S40_undefined']), 1)
        self.assertAlmostEqual(float(run['Precision@k_groups']), 2.)
        self.assertAlmostEqual(float(run['Recall@max_k']), 60 / 61)
        self.assertEqual(int(run['Recall@max_k_defined']), 1)


REPO = EXECUTION.parents[2]
FOUNDATION_API = REPO / 'stage2/01_bases/output/api'
sys.path.insert(0, str(FOUNDATION_API))
import hashlib  # noqa: E402
import pandas as pd  # noqa: E402
from bases.corpora import write_h5ad  # noqa: E402
from bases.fold_seed import build_fold_seed  # noqa: E402
from bases.significance import bh_adjust, group_key, membership, publish_statistics  # noqa: E402
import plot_tables  # noqa: E402
from evaluation_driver import run_driver  # noqa: E402

CASES4 = ['tenv:temp_env::Plate1-A01', 'tenv:temp_env::Plate1-A02',
          'tenv:temp_env::Plate1-A03', 'tenv:temp_env::Plate1-A04']
PARENTS = ['PARENT-A', 'PARENT-B', 'PARENT-C', 'PARENT-D']
FOLD_ROLES = {0: ['test', 'test', 'train', 'validation'],
              1: ['validation', 'train', 'test', 'test']}


def build_fixture(root):
    """A 4-well, 2-fold, 120-gene... no: full 18,164-gene synthetic release."""
    root = Path(root)
    axis = np.arange(GENES, dtype=np.float64)
    measured = np.stack([GENES - axis, axis + 1.0, np.ones(GENES), np.ones(GENES)])
    group_ids = [group_key(parent, 1.0) for parent in PARENTS]
    rows = []
    for fold, roles in FOLD_ROLES.items():
        for case, parent, group, role in zip(CASES4, PARENTS, group_ids, roles):
            rows.append({'sample_id': case, 'drug_parent': parent, 'dose_um': 1.0,
                         'group_id': group, 'role': role, 'fold': fold})
    cohorts = pd.DataFrame(rows)
    (root / 'fold').mkdir(parents=True, exist_ok=True)
    cohorts.to_csv(root / 'fold/cohorts.csv', index=False)
    digest = hashlib.sha256((root / 'fold/cohorts.csv').read_bytes()).hexdigest()
    build_fold_seed(root / 'fold/cohorts.csv', root / 'fold',
                    accepted_sha256=digest)

    genes = pd.DataFrame({'gene_index': range(GENES),
                          'hgnc_id': [f'HGNC:{index}' for index in range(GENES)]})
    samples = pd.DataFrame([{'group_id': group, 'sample_id': case,
                             'source_profile_id': 'native-' + case,
                             'condition': 'treated'}
                            for case, group in zip(CASES4, group_ids)])
    groups = pd.DataFrame([{'group_id': group, 'drug_parent': parent, 'dose_um': 1.0,
                            'status': 'eligible', 'reason': '',
                            'sample_ids_sha256': membership([case])}
                           for parent, group, case in zip(PARENTS, group_ids, CASES4)])
    pvalue = np.ones((4, GENES))
    pvalue[0, :60] = 1e-8
    pvalue[0, 110] = 1e-8
    pvalue[1, 0] = 1e-8
    padj = np.vstack([bh_adjust(row) for row in pvalue])
    log2fc = np.zeros((4, GENES))
    log2fc[0, :60] = 3.0
    log2fc[0, 110] = 3.0
    log2fc[1, 0] = 3.0
    publish_statistics(root / 'truth', genes, groups, samples,
                       {'pvalue': pvalue, 'padj': padj, 'log2_fold_change': log2fc},
                       {'fixture': True}, fixture=True)

    frame = pd.DataFrame({'expression_row': range(4), 'plate': 'P1'}, index=CASES4)
    write_h5ad(root / 'corpus.h5ad', frame,
               pd.DataFrame({'hgnc_id': [f'HGNC:{index}' for index in range(GENES)]}),
               measured.astype('float32'), {('P1',): np.zeros(GENES)}, 'tenv')
    return {'measured': measured, 'group_ids': group_ids}


def build_runs(root):
    """PRnet predicts the fixture truth; Mean predicts a constant zero Δ."""
    axis = np.arange(GENES, dtype=np.float64)
    descending = GENES - axis
    predictor = np.stack([descending, descending])        # A = truth, B = reversed order
    history = [{'epoch': 1, 'training_loss': 0.5, 'validation_loss': 0.4},
               {'epoch': 2, 'training_loss': 0.4, 'validation_loss': 0.3}]
    write_run(root, 'prnet', 0, 17, cases=[CASES[1], CASES[0]],
              prediction=predictor, history=history)
    write_run(root, 'mean', 0, 17, cases=[CASES[0], CASES[1]],
              prediction=np.zeros((2, GENES)))


class DriverEndToEnd(unittest.TestCase):
    def test_driver_scores_fold_writes_tables_and_both_plot_value_tables(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build_fixture(root)
            build_runs(root / 'runs')
            record = run_driver(root / 'runs', ['prnet', 'mean'], [0], [17],
                                root / 'fold', root / 'truth', root / 'corpus.h5ad',
                                root / 'output', root / 'figure', allow_fixture=True)
            for name in ('per_well.csv', 'per_group.csv', 'per_run.csv',
                         'metrics_values.csv', 'learning_curve_values.csv',
                         'evaluation.json'):
                self.assertTrue((root / 'output' / name).is_file(), name)
                self.assertTrue((root / 'figure' / name).is_file()
                                or name in ('per_well.csv', 'per_group.csv', 'per_run.csv',
                                            'evaluation.json'), name)
            self.assertEqual(record['schema'], 'stage2-tenv-evaluation-v1')
            self.assertEqual(record['truth']['gene_count'], GENES)
            self.assertEqual(record['plot_value_schema']['columns'],
                             list(plot_tables.PLOT_VALUE_COLUMNS))

            runs = pd.read_csv(root / 'output/per_run.csv', keep_default_na=False)
            self.assertEqual(len(runs), 2)
            prnet = runs[runs.entry.eq('prnet')].iloc[0]
            mean = runs[runs.entry.eq('mean')].iloc[0]
            self.assertAlmostEqual(float(prnet['PCC']), 0.)
            self.assertAlmostEqual(float(prnet['R2']), -1.)
            self.assertEqual(int(prnet['max_k@P80S40_defined']), 1)
            self.assertEqual(mean['PCC'], '')
            self.assertEqual(int(mean['PCC_defined']), 0)
            self.assertEqual(int(mean['PCC_undefined']), 2)

            # The per-run macro mean equals the equal-weight mean of per-well values.
            wells = pd.read_csv(root / 'output/per_well.csv', keep_default_na=False)
            prnet_pcc = wells[wells.entry.eq('prnet') & wells.metric.eq('PCC')]
            self.assertAlmostEqual(float(prnet['PCC']),
                                   prnet_pcc.value.astype(float).mean())
            prnet_rmse = wells[wells.entry.eq('prnet') & wells.metric.eq('RMSE')]
            self.assertAlmostEqual(float(prnet['RMSE']),
                                   prnet_rmse.value.astype(float).mean())
            mean_pcc = wells[wells.entry.eq('mean') & wells.metric.eq('PCC')]
            self.assertTrue(mean_pcc.status.eq('undefined').all())
            self.assertTrue(mean_pcc.reason.eq('constant_prediction').all())
            mean_rmse = wells[wells.entry.eq('mean') & wells.metric.eq('RMSE')]
            self.assertTrue(mean_rmse.status.eq('defined').all())

            groups = pd.read_csv(root / 'output/per_group.csv', keep_default_na=False)
            prnet_max = groups[groups.entry.eq('prnet') & groups.metric.eq('max_k@P80S40')]
            defined = prnet_max[prnet_max.status.eq('defined')]
            undefined = prnet_max[prnet_max.status.eq('undefined')]
            self.assertEqual(float(defined.value.iloc[0]), 75.)
            self.assertEqual(undefined.reason.tolist(), ['no_qualifying_k'])
            prnet_recall = groups[(groups.entry.eq('prnet'))
                                  & (groups.metric.eq('Recall@max_k'))
                                  & (groups.status.eq('defined'))]
            self.assertAlmostEqual(float(prnet_recall.value.iloc[0]), 60 / 61)

            metrics = pd.read_csv(root / 'output/metrics_values.csv', dtype=str,
                                  keep_default_na=False)
            self.assertEqual(tuple(metrics.columns), plot_tables.PLOT_VALUE_COLUMNS)
            self.assertEqual(metrics.iloc[0]['record'], 'figure')
            panels = metrics[metrics.record.eq('panel')].panel.tolist()
            self.assertEqual(panels, list(plot_tables.PANEL_ORDER))
            series = metrics[metrics.record.eq('series')]
            self.assertEqual(len(series), 2 * len(plot_tables.PANEL_ORDER))
            precision = metrics[(metrics.panel.eq('Precision@k'))
                               & (metrics.role.eq('curve_mean'))]
            first = precision[(precision.series.eq('prnet')) & (precision.x.eq('1.0'))]
            self.assertEqual(first.y.tolist(), ['1.0'])
            self.assertEqual(first.defined.tolist(), ['1'])
            mean_first = precision[(precision.series.eq('mean')) & (precision.x.eq('1.0'))]
            self.assertEqual(mean_first.y.tolist(), ['0.0'])
            spearman = metrics[(metrics.panel.eq('Spearman@k'))
                               & (metrics.role.eq('curve_mean'))]
            self.assertTrue(spearman[spearman.series.eq('mean')].defined.eq('0').all())
            self.assertTrue(spearman[spearman.series.eq('mean')].reason
                            .eq('no_defined_groups').all())
            bars = metrics[(metrics.panel.eq('PCC')) & (metrics.role.eq('bar_mean'))]
            self.assertEqual(len(bars), 2)
            self.assertEqual(bars[bars.series.eq('mean')].defined.tolist(), ['0'])

            curve = pd.read_csv(root / 'output/learning_curve_values.csv', dtype=str,
                                keep_default_na=False)
            self.assertEqual(tuple(curve.columns), plot_tables.PLOT_VALUE_COLUMNS)
            self.assertEqual(curve[curve.record.eq('panel')].panel.tolist(), ['prnet'])
            curve_series = curve[curve.record.eq('series')]
            self.assertEqual(sorted(curve_series.series),
                             ['fold_0/seed_17|train', 'fold_0/seed_17|validation'])
            self.assertEqual(curve[curve.record.eq('point')].shape[0], 4)
            self.assertEqual(record['learning_curve_skipped'], None)


def _rows_equal(left, right):
    if len(left) != len(right):
        return False
    for first, second in zip(left, right):
        for column in plot_tables.PLOT_VALUE_COLUMNS:
            a, b = first[column], second[column]
            if isinstance(a, float) and isinstance(b, float) and np.isnan(a) and np.isnan(b):
                continue
            if str(a) != str(b):
                return False
    return True


class SharedSchemaCompatibility(unittest.TestCase):
    """The #162 plotting API, once merged, must accept the driver's tables."""

    def test_fallback_matches_shared_plot_api_when_merged(self):
        shared = plot_tables.shared_plot_api()
        if shared is None:
            self.skipTest('shared ADR 0063 plotting API (#162) is not on this checkout')
        expected, measured, predicted, truth = synthetic_fold()
        score = score_run(expected, predicted, measured, truth, history=[
            {'epoch': 1, 'training_loss': 0.5, 'validation_loss': 0.4},
            {'epoch': 2, 'training_loss': 0.4, 'validation_loss': 0.3}])
        entries = figure_entries([score], ['run'])
        fallback = plot_tables._metrics_values_fallback(entries)
        shared_values = shared.metrics_values(entries)
        self.assertTrue(_rows_equal(plot_tables.plot_rows(fallback),
                                    shared.plot_rows(shared_values)))
        models = learning_models([score], ['run'])
        fallback_curve = plot_tables._learning_values_fallback(models)
        shared_curve = shared.learning_curve_values(models)
        self.assertTrue(_rows_equal(plot_tables.plot_rows(fallback_curve),
                                    shared.plot_rows(shared_curve)))


if __name__ == "__main__":
    unittest.main()
