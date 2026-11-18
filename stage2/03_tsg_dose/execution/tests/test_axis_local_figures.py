"""Visible data and zero in the CSV-backed Round 1 y-axis contract."""
from pathlib import Path
import sys
import json
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / 'figure/script'
sys.path.insert(0, str(SCRIPT))
sys.path.insert(0, str(SCRIPT.parents[1] / 'execution/axis_local'))
from render_axis_local import AxisScale


class AxisReference(unittest.TestCase):
    def test_broken_linear_axis_retains_uncertainty_and_zero(self):
        axis = AxisScale.choose([.61, .7, .8])
        self.assertEqual(len(axis.segments), 2)
        rendered = axis.map([0, .61, .7, .8])
        self.assertTrue(all(0 <= pos <= 1 for pos in rendered))
        self.assertTrue(all(a < b for a, b in zip(rendered, rendered[1:])))

    def test_mixed_sign_and_log_axes_preserve_their_references(self):
        axis = AxisScale.choose([-.2, -.1, .4])
        self.assertEqual(len(axis.segments), 1)
        self.assertTrue(0 < axis.map([0])[0] < 1)
        axis = AxisScale.choose([.01, .1, 1], logarithmic=True)
        self.assertEqual([segment['scale'] for segment in axis.segments], ['linear', 'log'])
        zero, low, mid, high = axis.map([0, .01, .1, 1])
        self.assertLess(zero, low)
        self.assertAlmostEqual(mid - low, high - mid)


class CsvOnlyRoundTrip(unittest.TestCase):
    def test_renderer_reproduces_all_figures_from_plot_csvs_without_reports(self):
        from hashlib import sha256
        from tempfile import TemporaryDirectory
        import csv

        from render_axis_local import (
            METHODS, METRICS, PNG_NAMES, SCENARIOS, SOURCES, prepare, render,
        )

        def save(path, table):
            with path.open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(table[0]))
                writer.writeheader()
                writer.writerows(table)

        from publish_axis_local_report import accept_candidate, digest, prepare_candidate

        with TemporaryDirectory() as root:
            root = Path(root)
            matrix = root / 'formal'
            reports = matrix / 'reports'
            for scenario in SCENARIOS:
                report = reports / scenario
                report.mkdir(parents=True)
                folds, spread, doses, learning = [], [], [], []
                for i, source in enumerate(SOURCES):
                    for j, method in enumerate(METHODS):
                        base = dict(source=source, scenario=scenario, method=method, seed=42)
                        values = dict(pcc=.64 + j * .01, spearman=.55, rmse=.03, r2=-.2)
                        for fold in range(3):
                            folds.append(dict(base, fold=fold, **values))
                            for epoch in range(1, 3):
                                learning.append(dict(base, fold=fold, epoch=epoch,
                                                     train_loss=5 + (1 + j / 50 + fold / 5) / epoch,
                                                     validation_mse=.1 / epoch + j / 500 + fold / 200))
                        for metric in METRICS:
                            spread.append(dict(base, metric=metric, mean=values[metric],
                                               sd=.005, n_defined_folds=3))
                        for dose in ((.05, .5, 5.) if i == 0 else (.2, 1., 10.)):
                            doses.append(dict(base, dose_um=dose, **values))
                save(report / 'folds.csv', folds)
                save(report / 'fold_spread.csv', spread)
                save(report / 'learning.csv', learning)
                if scenario == 'cold-drug':
                    save(report / 'per_dose.csv', doses)
            plot_csv = root / 'plot-csv'
            prepare(reports, plot_csv)
            seal_fixture(matrix, digest)
            receipt = prepare_candidate(matrix, root / 'first')
            self.assertTrue(accept_candidate(root / 'first', receipt, visually_inspected=True)['visual_inspection'])
            with (reports / 'cold-drug/folds.csv').open('a') as stream:
                stream.write('tampered after verification')
            with self.assertRaisesRegex(ValueError, 'accepted scientific artifact changed'):
                accept_candidate(root / 'first', receipt, visually_inspected=True)
            matrix.rename(root / 'hidden-matrix')
            render(plot_csv, root / 'second')
            self.assertEqual({path.name for path in (root / 'first').glob('*.png')}, set(PNG_NAMES))
            for name in PNG_NAMES:
                first = (root / 'first' / name).read_bytes()
                second = (root / 'second' / name).read_bytes()
                self.assertEqual(sha256(first).digest(), sha256(second).digest(), name)

    def test_report_candidate_requires_combined_scientific_acceptance(self):
        from tempfile import TemporaryDirectory

        from publish_axis_local_report import prepare_candidate

        with TemporaryDirectory() as root:
            root = Path(root)
            with self.assertRaisesRegex(ValueError, '84-origin scientific acceptance'):
                prepare_candidate(root / 'formal', root / 'candidate')
            self.assertFalse((root / 'candidate').exists())



def seal_fixture(matrix, digest):
    from render_axis_local import METHODS, SCENARIOS, SOURCES

    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return digest(path)

    plan_sha256 = write(matrix / 'formal_plan.json', {'fixture': True})
    scenario_hashes = {}
    for scenario in SCENARIOS:
        origins = []
        for source in SOURCES:
            for fold in range(3):
                selection = write(matrix / 'selection' / source / scenario / f'fold-{fold}.json', {})
                for method in METHODS:
                    origin = matrix / 'runs' / source / scenario / f'fold-{fold}' / method
                    manifest = write(origin / 'manifest.json', {})
                    origins.append(dict(source=source, scenario=scenario, fold=fold, method=method,
                                        seed=42, run_manifest_sha256=manifest, selection_sha256=selection))
        report = matrix / 'reports' / scenario
        scenario_hashes[scenario] = write(report / 'acceptance.json', dict(
            schema='tsg-dose-round1-axis-local-scenario-acceptance-v1', scenario=scenario,
            formal_fits=42, plan_sha256=plan_sha256, origins=origins,
            outputs={path.name: digest(path) for path in report.glob('*.csv')},
        ))
    write(matrix / 'formal_acceptance.json', dict(
        schema='tsg-dose-round1-axis-local-formal-acceptance-v1', formal_fits=84,
        plan_sha256=plan_sha256, scenario_receipt_sha256=scenario_hashes,
    ))
