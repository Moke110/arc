"""Visible data and zero in the CSV-backed Round 1 y-axis contract."""
from pathlib import Path
import sys
import json
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'axis_local'
sys.path.insert(0, str(SCRIPT))
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

    def test_log_curve_disjoins_zero_segment_without_dropping_zero(self):
        import math
        from matplotlib import pyplot as plt

        axis = AxisScale.choose([0, .1, .3, 0], logarithmic=True)
        fig, ax = plt.subplots()
        self.addCleanup(plt.close, fig)
        lines = axis.plot(ax, [1, 2, 3, 4], [0, .1, .3, 0])
        self.assertEqual(len(lines), 2)
        zero, positive = [line.get_ydata() for line in lines]
        self.assertTrue(math.isfinite(zero[0]) and math.isfinite(zero[3]))
        self.assertTrue(math.isnan(positive[0]) and math.isnan(positive[3]))
        self.assertTrue(all(math.isfinite(point) for point in positive[1:3]))


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
            accepted = accept_candidate(root / 'first', receipt, visually_inspected=True)
            self.assertTrue(accepted['visual_inspection'])
            self.assertEqual(accepted['resources']['ram_samples'], 1)
            self._assert_install_keeps_previous_release(root, accepted)
            with (reports / 'cold-drug/selected_pooled.csv').open('a') as stream:
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


    def _assert_install_keeps_previous_release(self, root, accepted):
        import shutil
        from install_axis_local_report import CHECKS, HEALTH_SCHEMA, digest, install, recover, safe_tree

        repo = root / 'test-checkout'
        step = repo / 'stage2/07_tsg_dose'
        active = step / 'figure'
        active.mkdir(parents=True)
        for number in range(19):
            (active / f'legacy-{number}.txt').write_text(f'preserved byte {number}')
        old = safe_tree(active)
        (step / 'execution/historical').mkdir(parents=True)
        (step / 'log').mkdir()
        seal = step / 'log/242-round1-report-acceptance.json'
        seal.write_text(json.dumps(dict(status='passed', files=old)))
        accepted_path = root / 'accepted.json'
        accepted_path.write_text(json.dumps(accepted))
        (repo / 'script').mkdir()
        health_path = root / 'health.json'
        checks = {}
        for name in CHECKS:
            script, log = repo / 'script' / name, root / f'{name}.log'
            script.write_text('echo pass')
            log.write_text('pass')
            checks[name] = dict(returncode=0, script_sha256=digest(script), log_sha256=digest(log))
        health_path.write_text(json.dumps(dict(schema=HEALTH_SCHEMA, commit='0' * 40, checks=checks)))
        candidate = root / 'release-candidate'
        shutil.copytree(root / 'first', candidate)
        evidence = install(step, candidate, accepted_path, health_path)
        self.assertEqual(evidence['previous_report_files'], old)
        self.assertEqual(safe_tree(step / 'execution/historical/report_126'), old)
        self.assertEqual(safe_tree(active), accepted['files'])
        self.assertFalse(candidate.exists())
        recover(step, candidate, accepted_path)
        self._assert_recover_verifies_successor(step, candidate, accepted_path)

    def _assert_recover_verifies_successor(self, step, candidate, accepted_path):
        import shutil
        from install_axis_local_report import AXIS_LOCAL_ARCHIVE, digest, recover

        source = SCRIPT.parents[1]
        successor = step / 'figure/baselines'
        shutil.copytree(source / AXIS_LOCAL_ARCHIVE / 'baselines', successor)
        with self.assertRaisesRegex(ValueError, 'successor receipt'):
            recover(step, candidate, accepted_path)

        receipt_path = step / 'log/round1-baselines-acceptance.json'
        original = json.loads((source / 'log/round1-baselines-acceptance.json').read_text())
        original['historical_receipts'] = {name: digest(step / 'log' / name)
                                          for name in original['historical_receipts']}
        receipt_path.write_text(json.dumps(original))
        recover(step, candidate, accepted_path)

        extra = step / 'figure/unsealed'
        extra.mkdir()
        with self.assertRaisesRegex(ValueError, 'installed report evidence changed'):
            recover(step, candidate, accepted_path)
        extra.rmdir()
        altered = successor / 'data/round1_pooled.csv'
        with altered.open('ab') as stream:
            stream.write(b'altered')
        with self.assertRaisesRegex(ValueError, 'successor report bytes'):
            recover(step, candidate, accepted_path)
        shutil.copyfile(source / AXIS_LOCAL_ARCHIVE / 'baselines/data/round1_pooled.csv', altered)

        receipt_path.write_text(json.dumps(dict(original, visual_inspection=False)))
        with self.assertRaisesRegex(ValueError, 'successor receipt'):
            recover(step, candidate, accepted_path)

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

    from publish_axis_local_report import TABLES, COLD_DRUG_TABLES
    from release_resources import PHASES

    from datetime import datetime, timezone
    stamp = datetime(2026, 11, 18, tzinfo=timezone.utc)
    (matrix / 'formal-cap-start.txt').write_text(
        'Id=arc-stage2-03-round1-formal.service\nInvocationID=' + 'a' * 32
        + '\nMemoryMax=10737418240\nMemorySwapMax=8589934592\nExecMainStartTimestamp=' + stamp.isoformat())
    (matrix / 'driver.time').write_text('Maximum resident set size (kbytes): 2000\nExit status: 0\n')
    for file in ('resources-10s.log', 'gpu-launch-window.log', 'gpu-10s.log'):
        (matrix / file).write_text(stamp.isoformat() + '\nMemTotal: 48000000 kB\n'
                                   'MemAvailable: 39000000 kB\nSwapFree: 64000000 kB\n2400, 10000, 10\n')
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
        tables = TABLES | (COLD_DRUG_TABLES if scenario == 'cold-drug' else set())
        for table in tables:
            if not (report / table).exists():
                (report / table).write_text('sealed ancillary table\n')
        (report / 'fit_resources.csv').write_text(
            'source,scenario,fold,method,seed,phase,action,started_utc_unix_seconds,seconds,returncode\n'
            + ''.join(f'{source},{scenario},{fold},{method},42,{phase},{action},{stamp.timestamp()},1,0\n'
                      for source in SOURCES for fold in range(3) for method in METHODS
                      for phase, action in PHASES))
        scenario_hashes[scenario] = write(report / 'acceptance.json', dict(
            schema='tsg-dose-round1-axis-local-scenario-acceptance-v1', scenario=scenario,
            formal_fits=42, plan_sha256=plan_sha256, origins=origins,
            outputs={path.name: digest(path) for path in report.glob('*.csv')},
        ))
    write(matrix / 'formal_acceptance.json', dict(
        schema='tsg-dose-round1-axis-local-formal-acceptance-v1', formal_fits=84,
        plan_sha256=plan_sha256, scenario_receipt_sha256=scenario_hashes,
    ))


class ReleaseRecovery(unittest.TestCase):
    def test_recovery_restores_sealed_trees_at_either_rename_boundary(self):
        from tempfile import TemporaryDirectory
        from install_axis_local_report import recover_swap, safe_tree

        for candidate_renamed in (False, True):
            with self.subTest(candidate_renamed=candidate_renamed), TemporaryDirectory() as temporary:
                root = Path(temporary)
                active, archive, candidate = (root / name for name in ('active', 'archive', 'candidate'))
                active.mkdir()
                (active / 'old.png').write_bytes(b'old report')
                candidate.mkdir()
                (candidate / 'new.png').write_bytes(b'new report')
                old, new = safe_tree(active), safe_tree(candidate)
                active.rename(archive)
                if candidate_renamed:
                    candidate.rename(active)
                recover_swap(active, archive, candidate, old, new)
                self.assertEqual(safe_tree(active), old)
                self.assertEqual(safe_tree(candidate), new)
                self.assertFalse(archive.exists())

    def test_telemetry_gap_cannot_pass_publication(self):
        from release_resources import verify_window

        with self.assertRaisesRegex(ValueError, '10-second GPU telemetry'):
            verify_window({0: {}, 10: {}, 40: {}}, 1, 39, 'GPU')


class SuccessorPresentation(unittest.TestCase):
    def test_legacy_inventory_requires_separate_receipt_for_tenv_kpgt(self):
        from shutil import copytree
        from tempfile import TemporaryDirectory

        from install_axis_local_report import (
            GRAY_ARCHIVE, digest, verify_baseline_presentation,
        )

        original = SCRIPT.parents[1]
        with TemporaryDirectory() as root:
            step = Path(root)
            copytree(original / 'log', step / 'log', ignore=lambda path, names: [
                name for name in names if name not in (
                    'baselines-gray-acceptance.json', 'round1-baselines-acceptance.json')
            ])
            active = step / 'figure'
            copytree(original / GRAY_ARCHIVE, active)
            verify_baseline_presentation(step, active)
            successor = active / 'tenv_kpgt'
            successor.mkdir()
            report = successor / 'metrics.png'
            report.write_bytes(b'accepted successor artifact')
            with self.assertRaisesRegex(ValueError, 'presentation bytes changed'):
                verify_baseline_presentation(step, active)
            receipt = step / 'log/275-kpgt-formal/acceptance.json'
            receipt.parent.mkdir()
            receipt.write_text(json.dumps(dict(
                schema='tsg-dose-tenv-kpgt-report-v1', status='passed', visual_inspection=True,
                files={'metrics.png': digest(report)},
            )))
            verify_baseline_presentation(step, active)
            report.write_bytes(b'changed successor artifact')
            with self.assertRaisesRegex(ValueError, 'presentation bytes changed'):
                verify_baseline_presentation(step, active)
            report.write_bytes(b'accepted successor artifact')
            (active / 'unsealed').mkdir()
            with self.assertRaisesRegex(ValueError, 'presentation bytes changed'):
                verify_baseline_presentation(step, active)

    def test_promoted_kpgt_and_archived_gray_match_their_receipts(self):
        from install_axis_local_report import verify_baseline_presentation

        step = SCRIPT.parents[1]
        verify_baseline_presentation(step, step / 'figure')
