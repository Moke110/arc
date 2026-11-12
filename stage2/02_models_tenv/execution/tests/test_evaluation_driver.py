"""Fail-closed run-set and @d scoring contract of the tenv evaluation driver (#168)."""
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

EXECUTION = Path(__file__).resolve().parents[1]
SCRIPT = EXECUTION / 'script'
sys.path.insert(0, str(SCRIPT))
from evaluation_driver import (  # noqa: E402
    CONFUSION_CLASSES, DIRECT_COLUMNS, DriverError, PER_GROUP_FIELDS,
    PER_RUN_FIELDS, confusion_rows, figure_entries, learning_models,
    load_run_set, per_group_rows, per_run_curve_rows, per_run_rows,
    per_well_rows, score_run, validate_run)

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
    """A known-value fold for the direct metrics: two wells, two groups."""
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


def reliability_fold():
    """A hand-computable two-group fold: one changed group, one genuine null.

    Group A predicts ``[4, 2, 1, 0, ...]`` and has three measured effects with
    ``|log2FC|`` ``[2, 1, 0.5]``. Group B predicts an all-zero Δ, so it is a null
    result. ``M = 4`` and the grid is ``d = 0.04 k``; at ``k = 50`` group A's
    window is ``{gene0, gene1}`` (two hits of three measured effects) and group B
    keeps an empty window, so the macro Reliability is ``(0.8 + 1) / 2``.
    """
    measured = np.zeros((2, GENES))
    measured[0, 0], measured[0, 1] = 2.0, 1.0
    measured[1, 0] = 5.0
    predicted = np.zeros((2, GENES))
    predicted[0, 0], predicted[0, 1], predicted[0, 2] = 4.0, 2.0, 1.0
    padj = np.full((2, GENES), 0.5)
    log2fc = np.zeros((2, GENES))
    padj[0, :3] = 1e-8
    log2fc[0, :3] = [2.0, 1.0, 0.5]
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


def row_at(rows, group_id, k):
    return next(row for row in rows if row['group_id'] == group_id and row['k'] == k)


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

    def test_per_group_curve_matches_hand_computed_reliability(self):
        expected, measured, predicted, truth = reliability_fold()
        score = score_run(expected, predicted, measured, truth)
        rows = per_group_rows([score])
        at_star_k = row_at(rows, 'group-A', 1)
        self.assertEqual(at_star_k['n_window'], 3)
        self.assertEqual(at_star_k['n_measured'], 3)
        self.assertEqual(at_star_k['n_hits'], 3)
        for metric in ('Precision@d', 'Recall@d', 'F1@d', 'Spearman@d', 'Reliability@d'):
            self.assertAlmostEqual(at_star_k[metric], 1.)
        middle = row_at(rows, 'group-A', 50)
        self.assertEqual(middle['n_window'], 2)
        self.assertEqual(middle['n_measured'], 3)
        self.assertEqual(middle['n_hits'], 2)
        self.assertAlmostEqual(middle['Precision@d'], 1.)
        self.assertAlmostEqual(middle['Recall@d'], 2 / 3)
        self.assertAlmostEqual(middle['F1@d'], 0.8)
        self.assertAlmostEqual(middle['Spearman@d'], 1.)
        self.assertAlmostEqual(middle['Reliability@d'], 0.8)
        strictest = row_at(rows, 'group-A', 100)
        self.assertEqual(strictest['n_window'], 1)
        self.assertAlmostEqual(strictest['F1@d'], 0.5)
        self.assertAlmostEqual(strictest['Reliability@d'], 0.5)

    def test_null_group_keeps_the_declared_totals_at_every_threshold(self):
        expected, measured, predicted, truth = reliability_fold()
        score = score_run(expected, predicted, measured, truth)
        rows = per_group_rows([score])
        null_rows = [row for row in rows if row['group_id'] == 'group-B']
        self.assertEqual(len(null_rows), 100)
        for row in null_rows:
            self.assertEqual(row['n_window'], 0)
            self.assertEqual(row['n_measured'], 0)
            self.assertEqual(row['n_hits'], 0)
            for metric in ('Precision@d', 'Recall@d', 'F1@d', 'Spearman@d', 'Reliability@d'):
                self.assertAlmostEqual(row[metric], 1.)
        run = per_run_rows([score])[0]
        self.assertEqual(run['groups'], 2)
        self.assertEqual(run['null_groups'], 1)
        self.assertEqual(run['changed_groups'], 1)
        self.assertAlmostEqual(run['null_Reliability@d*'], 1.)

    def test_d_star_is_the_argmax_of_the_recorded_macro_curve(self):
        expected, measured, predicted, truth = reliability_fold()
        score = score_run(expected, predicted, measured, truth)
        run = per_run_rows([score])[0]
        curve = per_run_curve_rows([score])
        reliability = [row['Reliability@d'] for row in curve]
        self.assertEqual(run['d_star_k'], int(np.argmax(reliability)) + 1)
        best = curve[run['d_star_k'] - 1]
        self.assertAlmostEqual(run['Reliability@d*'], best['Reliability@d'])
        self.assertAlmostEqual(run['Precision@d*'], best['Precision@d'])
        self.assertAlmostEqual(run['Recall@d*'], best['Recall@d'])
        self.assertAlmostEqual(run['F1@d*'], best['F1@d'])
        self.assertAlmostEqual(run['Spearman@d*'], best['Spearman@d'])
        self.assertEqual(run['d_star_k'], 1)
        self.assertAlmostEqual(run['d_star'], 0.04)
        self.assertAlmostEqual(run['d_star_prime'], 0.01)
        self.assertAlmostEqual(run['M'], 4.0)
        self.assertAlmostEqual(curve[49]['Reliability@d'], 0.9)
        self.assertAlmostEqual(curve[49]['null_Reliability@d'], 1.0)
        self.assertAlmostEqual(curve[49]['changed_Reliability@d'], 0.8)

    def test_per_group_table_covers_every_group_at_every_threshold(self):
        expected, measured, predicted, truth = reliability_fold()
        score = score_run(expected, predicted, measured, truth)
        rows = per_group_rows([score])
        self.assertEqual(len(rows), 2 * 100)
        for group_id in ('group-A', 'group-B'):
            group_rows = [row for row in rows if row['group_id'] == group_id]
            self.assertEqual([row['k'] for row in group_rows], list(range(1, 101)))
            for row in group_rows:
                self.assertAlmostEqual(row['d_prime'], row['k'] / 100)
                self.assertAlmostEqual(row['d'], row['d_prime'] * 4.0)
                self.assertEqual(row['n_wells'], 1)
        self.assertEqual({row['drug_parent'] for row in rows}, {'PARENT-A', 'PARENT-B'})

    def test_confusion_counts_classify_every_gene_at_d_star(self):
        expected, measured, predicted, truth = reliability_fold()
        score = score_run(expected, predicted, measured, truth)
        rows = confusion_rows([score])
        self.assertEqual(len(rows), len(CONFUSION_CLASSES) ** 2)
        counts = {(row['predicted_class'], row['measured_class']): row['count']
                  for row in rows}
        self.assertEqual(counts[('up', 'up')], 3)
        self.assertEqual(counts[('not_significant', 'not_significant')], 2 * GENES - 3)
        self.assertEqual(sum(counts.values()), 2 * GENES)
        self.assertEqual(rows[0]['d_star_k'], 1)

    def test_run_tables_record_the_grid_identity_and_direct_metrics(self):
        expected, measured, predicted, truth = synthetic_fold()
        score = score_run(expected, predicted, measured, truth)
        run = per_run_rows([score])[0]
        self.assertEqual(run['wells'], 2)
        self.assertEqual(run['threshold_points'], 100)
        self.assertEqual(run['threshold_rule'], 'd = k/100 * M')
        self.assertEqual(tuple(per_run_rows([score])[0]), tuple(PER_RUN_FIELDS))
        self.assertAlmostEqual(run['PCC'], 0.)
        self.assertAlmostEqual(run['R2'], -1.)
        self.assertEqual(run['PCC_defined_wells'], 2)
        self.assertEqual(run['PCC_undefined_wells'], 0)

    def test_zero_predicted_magnitude_cannot_define_the_scan(self):
        expected, measured, _, truth = reliability_fold()
        with self.assertRaisesRegex(DriverError, 'threshold scan'):
            score_run(expected, np.zeros((2, GENES)), measured, truth)

    def test_non_eligible_truth_raises_instead_of_dropping_a_group(self):
        expected, measured, predicted, truth = reliability_fold()
        truth.rows['group-B']['status'] = 'ineligible'
        with self.assertRaisesRegex(DriverError, 'not eligible'):
            score_run(expected, predicted, measured, truth)


REPO = EXECUTION.parents[2]
FOUNDATION_API = REPO / 'stage2/01_bases/output/api'
sys.path.insert(0, str(FOUNDATION_API))
from bases.corpora import write_h5ad  # noqa: E402
from bases.fold_seed import build_fold_seed  # noqa: E402
from bases.significance import bh_adjust, group_key, membership, publish_statistics  # noqa: E402
import plot_tables  # noqa: E402
from evaluation_driver import run_driver, verify_identity_hashes  # noqa: E402

CASES4 = ['tenv:temp_env::Plate1-A01', 'tenv:temp_env::Plate1-A02',
          'tenv:temp_env::Plate1-A03', 'tenv:temp_env::Plate1-A04']
PARENTS = ['PARENT-A', 'PARENT-B', 'PARENT-C', 'PARENT-D']
FOLD_ROLES = {0: ['test', 'test', 'train', 'validation'],
              1: ['validation', 'train', 'test', 'test']}


def build_fixture(root):
    """A 4-well, 2-fold, full 18,164-gene synthetic release."""
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
    """PRnet predicts the fixture truth; Mean predicts a constant nonzero Δ."""
    axis = np.arange(GENES, dtype=np.float64)
    descending = GENES - axis
    predictor = np.stack([descending, descending])        # A = truth, B = reversed order
    history = [{'epoch': 1, 'training_loss': 0.5, 'validation_loss': 0.4},
               {'epoch': 2, 'training_loss': 0.4, 'validation_loss': 0.3}]
    write_run(root, 'prnet', 0, 17, cases=[CASES[1], CASES[0]],
              prediction=predictor, history=history)
    write_run(root, 'mean', 0, 17, cases=[CASES[0], CASES[1]],
              prediction=np.full((2, GENES), 0.5))


def plot_values(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


class DriverEndToEnd(unittest.TestCase):
    def driver(self, root):
        record = run_driver(root / 'runs', ['prnet', 'mean'], [0], [17],
                            root / 'fold', root / 'truth', root / 'corpus.h5ad',
                            root / 'output', root / 'figure', allow_fixture=True)
        return record

    def test_driver_writes_the_new_tables_and_plot_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build_fixture(root)
            build_runs(root / 'runs')
            record = self.driver(root)
            for name in ('per_well.csv', 'per_group.csv', 'per_run.csv',
                         'per_run_curve.csv', 'confusion.csv', 'metrics_values.csv',
                         'metrics_indirect_values.csv', 'learning_curve_values.csv',
                         'evaluation.json'):
                self.assertTrue((root / 'output' / name).is_file(), name)
            for name in ('metrics_values.csv', 'metrics_indirect_values.csv',
                         'learning_curve_values.csv'):
                self.assertTrue((root / 'figure' / name).is_file(), name)
            self.assertEqual(record['schema'], 'stage2-tenv-evaluation-v2')
            self.assertEqual(record['truth']['gene_count'], GENES)
            self.assertEqual(record['threshold_grid']['points'], 100)
            self.assertEqual(record['plot_value_schema']['columns'],
                             list(plot_tables.PLOT_VALUE_COLUMNS))
            self.assertEqual(record['direct_metrics'], list(DIRECT_COLUMNS))
            self.assertEqual(record['reliability_metrics'],
                             ['Precision@d', 'Recall@d', 'F1@d', 'Spearman@d',
                              'Reliability@d'])

    def test_tables_cover_every_test_group_and_use_only_the_new_schema(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build_fixture(root)
            build_runs(root / 'runs')
            self.driver(root)
            groups = pd.read_csv(root / 'output/per_group.csv', keep_default_na=False)
            self.assertEqual(tuple(groups.columns), tuple(PER_GROUP_FIELDS))
            self.assertEqual(len(groups), 2 * 2 * 100)
            self.assertFalse([column for column in groups.columns
                              if '@k' in column or 'AUPRC' in column])
            runs = pd.read_csv(root / 'output/per_run.csv', keep_default_na=False)
            self.assertEqual(tuple(runs.columns), tuple(PER_RUN_FIELDS))
            self.assertEqual(len(runs), 2)
            curve = pd.read_csv(root / 'output/per_run_curve.csv', keep_default_na=False)
            self.assertEqual(len(curve), 2 * 100)
            confusion = pd.read_csv(root / 'output/confusion.csv', keep_default_na=False)
            self.assertEqual(len(confusion), 2 * 9)
            self.assertEqual(int(confusion['count'].astype(int).sum()), 2 * 2 * GENES)
            record = json.loads((root / 'output/evaluation.json').read_text())
            self.assertNotIn('k_axis', record)
            self.assertNotIn('indirect_metrics', record)
            for run in record['runs']:
                self.assertNotIn('n_gated_truth', run)
                self.assertGreater(float(run['M']), 0)
                self.assertEqual(run['threshold_points'], 100)

    def test_run_summary_d_star_agrees_with_the_recorded_macro_curve(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build_fixture(root)
            build_runs(root / 'runs')
            record = self.driver(root)
            runs = pd.read_csv(root / 'output/per_run.csv', keep_default_na=False)
            curve = pd.read_csv(root / 'output/per_run_curve.csv', keep_default_na=False)
            for record_run in record['runs']:
                key = (record_run['entry'], record_run['fold'], record_run['seed'])
                row = runs[(runs.entry == key[0]) & (runs.fold.astype(int) == key[1])
                           & (runs.seed.astype(int) == key[2])].iloc[0]
                own = curve[(curve.entry == key[0]) & (curve.fold.astype(int) == key[1])
                            & (curve.seed.astype(int) == key[2])]
                reliability = own['Reliability@d'].astype(float).to_numpy()
                self.assertEqual(int(row['d_star_k']), int(np.argmax(reliability)) + 1)
                best = own.iloc[int(row['d_star_k']) - 1]
                self.assertAlmostEqual(float(row['Reliability@d*']),
                                       float(best['Reliability@d']))
                self.assertAlmostEqual(float(row['M']), float(record_run['M']))
                self.assertAlmostEqual(float(row['d_star']), float(best['d']))
                self.assertAlmostEqual(float(row['d_star']),
                                       int(row['d_star_k']) / 100 * float(row['M']))
                self.assertEqual(record_run['d_star_k'], int(row['d_star_k']))
                self.assertAlmostEqual(record_run['Reliability@d*'],
                                       float(row['Reliability@d*']))

    def test_run_macro_equals_the_mean_of_its_group_curves(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build_fixture(root)
            build_runs(root / 'runs')
            self.driver(root)
            groups = pd.read_csv(root / 'output/per_group.csv', keep_default_na=False)
            curve = pd.read_csv(root / 'output/per_run_curve.csv', keep_default_na=False)
            for (entry, fold, seed), own in curve.groupby(['entry', 'fold', 'seed']):
                members = groups[groups.entry.eq(entry)
                                 & groups.fold.astype(int).eq(int(fold))
                                 & groups.seed.astype(int).eq(int(seed))]
                for k in (1, 50, 100):
                    macro = float(own[own.k.astype(int).eq(k)]['Reliability@d'].iloc[0])
                    group_mean = float(
                        members[members.k.astype(int).eq(k)]['Reliability@d']
                        .astype(float).mean())
                    self.assertAlmostEqual(macro, group_mean, places=12)
                null_groups = members[members.n_measured.astype(int).eq(0)].group_id.nunique()
                null_rows = members[members.n_measured.astype(int).eq(0)]
                self.assertEqual(len(null_rows), 100 * null_groups)

    def test_scoring_uses_only_the_evaluated_folds_test_groups(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = build_fixture(root)
            build_runs(root / 'runs')
            self.driver(root)
            groups = pd.read_csv(root / 'output/per_group.csv', keep_default_na=False)
            self.assertEqual(set(groups.group_id), set(fixture['group_ids'][:2]))
            self.assertEqual(groups.fold.astype(int).unique().tolist(), [0])

    def test_direct_and_reliability_plot_values_carry_the_new_metric_set(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build_fixture(root)
            build_runs(root / 'runs')
            self.driver(root)
            metrics = plot_values(root / 'output/metrics_values.csv')
            self.assertEqual(tuple(metrics.columns), plot_tables.PLOT_VALUE_COLUMNS)
            panels = metrics[metrics.record.eq('panel')].panel.tolist()
            self.assertEqual(panels, list(plot_tables.DIRECT_PANEL_ORDER))
            self.assertEqual(metrics.iloc[0]['kind'], 'metrics')
            indirect = plot_values(root / 'output/metrics_indirect_values.csv')
            panels = indirect[indirect.record.eq('panel')].panel.tolist()
            self.assertEqual(panels, list(plot_tables.INDIRECT_PANEL_ORDER))
            self.assertIn('Reliability@d', panels)
            self.assertIn('Confusion@d*', panels)
            self.assertFalse(indirect[['panel', 'series']].apply(
                lambda column: column.str.contains('@k', regex=False)).any().any())
            markers = indirect[indirect.role.eq('marker')]
            self.assertEqual(len(markers), 2)
            reference = indirect[indirect.role.eq('reference')]
            self.assertEqual(len(reference), 2 * 100)
            heat = indirect[indirect.role.eq('heat_cell')]
            self.assertEqual(len(heat), 2 * 9)
            self.assertEqual(int(heat['y'].astype(float).sum()), 2 * 2 * GENES)

            curve = plot_values(root / 'output/learning_curve_values.csv')
            self.assertEqual(curve[curve.record.eq('panel')].panel.tolist(), ['prnet'])
            self.assertEqual(sorted(curve[curve.record.eq('series')].series),
                             ['fold_0/seed_17|train', 'fold_0/seed_17|validation'])

    def test_mean_run_keeps_constant_prediction_reasons(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build_fixture(root)
            build_runs(root / 'runs')
            self.driver(root)
            wells = pd.read_csv(root / 'output/per_well.csv', keep_default_na=False)
            mean_pcc = wells[wells.entry.eq('mean') & wells.metric.eq('PCC')]
            self.assertTrue(mean_pcc.status.eq('undefined').all())
            self.assertTrue(mean_pcc.reason.eq('constant_prediction').all())
            mean_rmse = wells[wells.entry.eq('mean') & wells.metric.eq('RMSE')]
            self.assertTrue(mean_rmse.status.eq('defined').all())
            # The per-run direct macro equals the equal-weight mean of per-well values.
            runs = pd.read_csv(root / 'output/per_run.csv', keep_default_na=False)
            prnet = runs[runs.entry.eq('prnet')].iloc[0]
            prnet_rmse = wells[wells.entry.eq('prnet') & wells.metric.eq('RMSE')]
            self.assertAlmostEqual(float(prnet['RMSE']),
                                   prnet_rmse.value.astype(float).mean())
            self.assertEqual(prnet['PCC_defined_wells'], 2)


class PlotValueSeam(unittest.TestCase):
    def test_indirect_values_require_one_shared_d_prime_axis(self):
        expected, measured, predicted, truth = reliability_fold()
        first = score_run(expected, predicted, measured, truth, entry='first')
        second = score_run(expected, predicted, measured, truth, entry='second')
        entries = figure_entries([first, second], ['first', 'second'])
        entries[1]['runs'][0]['d_prime'] = [0.5] * 100
        with self.assertRaisesRegex(ValueError, 'd_prime axis'):
            plot_tables.indirect_metrics_values(entries)

    def test_figure_entries_carry_the_threshold_marker_values(self):
        expected, measured, predicted, truth = reliability_fold()
        score = score_run(expected, predicted, measured, truth, history=[
            {'epoch': 1, 'training_loss': 0.5, 'validation_loss': 0.4}])
        entries = figure_entries([score], ['run'])
        run = entries[0]['runs'][0]
        self.assertEqual(run['d_star_k'], 1)
        self.assertAlmostEqual(run['reliability_at_d_star'], 1.0)
        self.assertEqual(run['confusion']['classes'], list(CONFUSION_CLASSES))
        models = learning_models([score], ['run'])
        self.assertEqual(models[0]['name'], 'run')


class IdentityDigestBinding(unittest.TestCase):
    """A completed run's bytes must match the digests its identity claims."""

    def _run(self, root, *, prediction_sha256='missing', case_sha256='missing'):
        run = root / 'prnet' / 'fold_0' / 'seed_17'
        run.mkdir(parents=True)
        np.save(run / 'prediction.npy', np.arange(4, dtype='float32'))
        (run / 'case_ids.json').write_text('["case"]')
        identity = {}
        if prediction_sha256 != 'missing':
            identity['prediction_sha256'] = prediction_sha256
        if case_sha256 != 'missing':
            identity['case_ids_sha256'] = case_sha256
        return SimpleNamespace(run_dir=run, identity=identity), run

    def test_matching_digests_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            run, _ = self._run(Path(temporary))
            run.identity['prediction_sha256'] = hashlib.sha256(
                (run.run_dir / 'prediction.npy').read_bytes()).hexdigest()
            run.identity['case_ids_sha256'] = hashlib.sha256(
                (run.run_dir / 'case_ids.json').read_bytes()).hexdigest()
            verify_identity_hashes(run)

    def test_missing_digest_fails_a_formal_run(self):
        with tempfile.TemporaryDirectory() as temporary:
            run, _ = self._run(Path(temporary))
            with self.assertRaisesRegex(DriverError, 'missing prediction_sha256'):
                verify_identity_hashes(run)

    def test_mismatched_digest_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            run, _ = self._run(Path(temporary), prediction_sha256='0' * 64)
            with self.assertRaisesRegex(DriverError, 'prediction digest differs'):
                verify_identity_hashes(run)

    def test_fixture_runs_may_omit_digests(self):
        with tempfile.TemporaryDirectory() as temporary:
            run, _ = self._run(Path(temporary))
            verify_identity_hashes(run, allow_fixture=True)


if __name__ == '__main__':
    unittest.main()
