"""Compact checks for selection, replay and budget disclosure without real inputs."""
from pathlib import Path
import os
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from xml.etree import ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from artifacts import read_json, rows, run_directory, sha256, write_csv, write_json
from deliver import publish
from learning_curves import V2_MODELS, render
from supplement_overlay import (ATTEMPT, EVALUATION_FILES, MODELS, SLOTS, TARGET,
                                evaluate, evaluation_proof, files, materialize,
                                pinned_json, selection, slot_key, verify_copy, verify_seal)
from test_delivery import evaluation_fixture


def v2_history(slot, iteration):
    return {'identity': dict(zip(('model', 'fold', 'seed'), slot)), 'protocol': V2_MODELS[slot[0]],
            'history': [dict(iteration=iteration, training_loss=2, validation_loss=1,
                             validation_conditions_evaluated=138)]}


def fixture(root, status='completed', original_history=True):
    original = []
    supplemental = root / 'supplement-cmonge-v2-fold0-seed17'
    source = root / 'formal'
    for slot in SLOTS:
        model, fold, seed = slot
        directory = run_directory(source, *slot)
        identity = dict(model=model, fold=fold, seed=seed)
        run_status = 'timeout' if slot == TARGET else 'completed'
        write_json(directory / 'run.json', dict(status=run_status, identity=identity))
        if slot == TARGET and original_history:
            write_json(directory / 'history.json', v2_history(slot, 50))
        original.append(dict(identity, identity=identity, status=run_status, command=None, directory=str(directory),
                             files=files(directory)))
    target = next(row for row in original if slot_key(row) == TARGET)
    target['budget'] = {'remaining_at_launch_seconds': 424.49, 'prior_debit_seconds': 13975.51,
                        'attempt_wall_seconds': 424.61, 'peak_host_rss_gib': None}
    directory = run_directory(supplemental, *TARGET)
    write_json(directory / 'run.json', dict(status=status, identity=v2_history(TARGET, 150)['identity']))
    if status == 'completed':
        write_json(directory / 'history.json', v2_history(TARGET, 150))
    selected = dict(target, directory=str(directory), status=status, files=files(directory))
    terminal = dict(status=status, attempt_id=ATTEMPT, run_sha256=selected['files']['run.json'],
                    account=dict(timing=dict(pre_inference_watchdog_seconds=800, adapter_training_seconds=770,
                                            observed_inference_seconds=40, adapter_inference_seconds=39),
                                 peak_host_rss_gib=2.2, time_peak_rss_gib=2.1))
    return dict(original_output=str(source), input_sha256='fixture', input_root=str(root / 'input'),
                original=original, terminal=terminal,
                selected=[selected if row is target else row for row in original])


def seal_fixture(output, record):
    entries = [{k: row[k] for k in ('model', 'fold', 'seed', 'status', 'command')}
               | {'run_sha256': row['files']['run.json']} for row in record['selected']]
    write_json(output / 'predictions.json', dict(input_sha256='fixture', models=list(MODELS),
               runs=entries, seeds=[17, 42, 2026], folds=[0, 1, 2], sealed=True))


class SupplementOverlayTests(unittest.TestCase):
    def test_terminal_authentication_binds_slot_status_and_original_timeout(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            record = fixture(root, 'timeout')
            formal = root / 'execution/formal'
            formal.mkdir(parents=True)
            write_json(root / 'input/manifest.json', {})
            original = next(r for r in record['original'] if slot_key(r) == TARGET)
            original.update(scheduler_sha256='scheduler', account_files={'account.json': 'account'})
            selected = next(r for r in record['selected'] if slot_key(r) == TARGET)
            timeout = dict(run_files=original['files'], scheduler_sha256='scheduler', summary_sha256='summary',
                           account_files=original['account_files'])
            directory = Path(selected['directory'])
            policy = dict(output_root=str(directory.parents[3]), original=timeout)
            terminal = dict(record['terminal'], identity=selected['identity'], run_directory=str(directory),
                            files=selected['files'], original=timeout, policy_sha256='policy')
            queue = dict(output_root=record['original_output'], input_root=record['input_root'])

            def select():
                return selection(formal, 'stage', root / 'policy.json', 'policy', 'terminal')

            def verify_terminal(*args):
                self.assertEqual(args[-1], 'terminal')
                return policy, terminal

            with (patch.dict(sys.modules, cmonge_supplement=SimpleNamespace(verify_terminal=verify_terminal)),
                  patch('supplement_overlay.preflight', return_value=(queue, {})),
                  patch('supplement_overlay.original_records', return_value=(record['original'], 'summary'))):
                outcome = select()
                self.assertEqual(len(outcome['selected']), 63)
                self.assertEqual([r.get('attempt_id') for r in outcome['selected']].count(ATTEMPT), 1)
                for field, invalid in [('status', 'running'), ('status', 'completed'), ('attempt_id', 'wrong'),
                                       ('identity', dict(terminal['identity'], seed=42)),
                                       ('run_directory', original['directory']), ('run_sha256', 'wrong')]:
                    with self.subTest(field=field), patch.dict(terminal, {field: invalid}):
                        with self.assertRaisesRegex(ValueError, 'terminal selection mismatch'):
                            select()
                with patch.dict(timeout, scheduler_sha256='forged'):
                    with self.assertRaisesRegex(ValueError, 'original timeout differs'):
                        select()

    def test_two_terminal_paths_select_exactly_one_slot_without_touching_source(self):
        for status in ('completed', 'timeout'):
            with self.subTest(status=status), tempfile.TemporaryDirectory() as temporary:
                record = fixture(Path(temporary), status)
                output = materialize(record)
                self.assertEqual(pinned_json(output / 'selection.json', sha256(output / 'selection.json')), record)
                verify_copy(output, record)
                actual = [slot_key(row) for row in record['selected']]
                self.assertEqual(actual, list(SLOTS))
                original = next(row for row in record['original'] if slot_key(row) == TARGET)
                selected = run_directory(output, *TARGET)
                self.assertEqual(read_json(selected / 'run.json')['status'], status)
                self.assertEqual(read_json(Path(original['directory']) / 'run.json')['status'], 'timeout')
                self.assertEqual(files(original['directory']), original['files'])
                source = next(Path(row['directory']) for row in record['selected'] if slot_key(row) == TARGET)
                self.assertNotEqual((selected / 'run.json').stat().st_ino, (source / 'run.json').stat().st_ino)
                with self.assertRaisesRegex(ValueError, 'already exists'):
                    materialize(record)

    def test_missing_extra_and_linked_files_fail_verification(self):
        with tempfile.TemporaryDirectory() as temporary:
            record = fixture(Path(temporary))
            output = materialize(record)
            target = run_directory(output, *TARGET)
            extra = target / 'extra'
            extra.write_text('extra')
            with self.assertRaisesRegex(ValueError, 'slot/file mismatch'):
                verify_copy(output, record)
            extra.unlink()
            original = Path(next(row['directory'] for row in record['selected'] if slot_key(row) == TARGET))
            path = target / 'history.json'
            path.unlink()
            with self.assertRaisesRegex(ValueError, 'slot/file mismatch'):
                verify_copy(output, record)
            path.symlink_to(original / path.name)
            with self.assertRaisesRegex(ValueError, 'linked path'):
                verify_copy(output, record)
            path.unlink()
            os.link(original / path.name, path)
            with self.assertRaisesRegex(ValueError, 'regular file'):
                verify_copy(output, record)

    def test_forged_seal_stops_before_frozen_evaluator(self):
        with tempfile.TemporaryDirectory() as temporary:
            record = fixture(Path(temporary), 'timeout')
            output = materialize(record)

            def forged_seal(*_):
                seal_fixture(output, record)
                path = output / 'predictions.json'
                seal = read_json(path)
                next(row for row in seal['runs'] if slot_key(row) == TARGET)['status'] = 'completed'
                write_json(path, seal)

            with patch('benchmark.freeze_predictions', side_effect=forged_seal), patch('benchmark.evaluate') as score:
                with self.assertRaisesRegex(ValueError, 'frozen seal differs'):
                    evaluate(output, record, sha256(output / 'selection.json'))
                score.assert_not_called()
            seal_fixture(output, record)
            verify_seal(output, record)
            with self.assertRaisesRegex(ValueError, 'reviewed digest'):
                pinned_json(output / 'selection.json', 'wrong')

    def test_publish_preserves_base_evidence_and_separates_attempt_histories(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            record = fixture(root)
            output = materialize(record)
            tables = root / 'tables'
            evaluation_fixture(tables)
            # A small artifact-only replay fixture with a completed supplemental row.
            for name in EVALUATION_FILES:
                (output / name).write_bytes((tables / name).read_bytes())
            for name in ('summary.csv', 'per_run.csv', 'pr_bands.csv', 'best_runs.csv'):
                table = rows(output / name)
                for row in table:
                    if row['model'] == 'failed_model':
                        row['model'] = 'cmonge'
                if name == 'per_run.csv':
                    row = next(r for r in table if slot_key(dict(r, fold=int(r['fold']), seed=int(r['seed']))) == TARGET)
                    row['status'] = 'completed'
                write_csv(output / name, table)
            seal_fixture(output, record)
            base = {name: sha256(output / name) for name in EVALUATION_FILES}
            write_json(output / 'evaluation.json', {'artifacts': base, 'workers': 2})
            selection_sha = sha256(output / 'selection.json')
            proof = evaluation_proof(output, record, selection_sha)
            destination = root / 'presentation'
            manifest = publish(output, destination, supplement=(record, proof, 'pinned-test'))
            for name, value in base.items():
                self.assertEqual(sha256(output / name), value)
            for name in ('metrics', 'pr_panels', 'best_run', 'learning_curves'):
                root_svg = ET.parse(destination / f'{name}.svg')
                self.assertIn('14,400 s', ' '.join(root_svg.getroot().itertext()))
            data = read_json(destination / 'learning_curves.json')
            self.assertEqual(data['terminal_outcomes'], 63)
            self.assertEqual(data['original_v2_records'][0]['validation'], [[50, 1]])
            selected = next(row for row in data['records'] if slot_key(row) == TARGET)
            self.assertEqual(selected['validation'], [[150, 1]])
            self.assertEqual(selected['attempt_id'], ATTEMPT)
            target = next(row for row in rows(destination / 'per_run.csv') if
                          (row['model'], row['fold'], row['seed']) == ('cmonge', '0', '17'))
            self.assertEqual(target['selected_attempt_id'], ATTEMPT)
            self.assertEqual(target['original_attempt_status'], 'timeout')
            self.assertEqual(manifest['source_artifacts'], base)
            self.assertNotEqual(manifest['artifacts']['per_run.csv'], base['per_run.csv'])
            with self.assertRaisesRegex(ValueError, 'authenticated supplement'):
                publish(output, destination)

    def test_missing_original_and_supplement_histories_are_reported_as_timeouts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            record = fixture(root, 'timeout', original_history=False)
            output = materialize(record)
            proof = dict(original_timeout=next(r for r in record['original'] if slot_key(r) == TARGET),
                         selected_attempt=record['terminal'], selection_sha256='fixture')
            data = render(output, root / 'presentation', final=True, supplement=(record, proof, 'fixture'))
            original, = data['original_v2_records']
            selected = next(row for row in data['records'] if slot_key(row) == TARGET)
            for attempt in (original, selected):
                self.assertEqual(attempt['status'], 'timeout')
                self.assertEqual(attempt['validation'], [])
                self.assertEqual(attempt['history_status'], 'missing')
            self.assertEqual(len(data['records']), 63)


if __name__ == '__main__':
    unittest.main()
