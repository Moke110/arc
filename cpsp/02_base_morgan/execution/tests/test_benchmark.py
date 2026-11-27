"""Complete synthetic matrix through the public command dispatcher, no formal fits."""
import contextlib
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import benchmark
import training
import workflow


class BenchmarkAcceptance(unittest.TestCase):
    def invoke(self, *args):
        with patch.object(sys, 'argv', ['workflow.py', *map(str, args)]), contextlib.redirect_stdout(io.StringIO()):
            workflow.main()

    def test_complete_60_fit_smoke_reuse_regeneration_and_corruption(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = [dict(condition=condition, phenotype='cell_viability', canonical_smiles=smiles,
                         dose_um=dose, direction=parent % 3 - 1, source_label=f'{parent}-{alias}',
                         source_asset='synthetic-fixture', source_cell=f'{parent}-{dose}-{alias}')
                    for condition in ['baseline', 'H-R']
                    for parent, smiles in enumerate(['C', 'CC', 'CCC', 'CCCC', 'CCO', 'CCN', 'CO', 'CN', 'CCCl'])
                    for dose in [0.04, 0.2, 1, 5] for alias in range(2)]
            source = root/'fixture.csv'
            with source.open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0])
                writer.writeheader()
                writer.writerows(rows)
            self.invoke('prepare', '--input', source, '--output', root/'prepared', '--settings', *workflow.SETTINGS)
            prepared = root/'prepared/prepared.json'
            output = root/'benchmark'
            calls = []
            def inline_guard(command, log, **kwargs):
                # Test-only process launcher replacement: real tiny models, no host telemetry claim.
                calls.append(command)
                args = command[command.index('train')+1:]
                args = dict(zip(args[::2], args[1::2]))
                training.train(args['--prepared'], args['--output'], args['--condition'], args['--phenotype'],
                               args['--setting'], int(args['--fold']), int(args['--seed']), int(args['--smoke-epochs']))
                (log/'time.txt').write_text('Maximum resident set size (kbytes): 1\n')
                return dict(exit_code=0, stop_reason=None, peak_tree_rss_bytes=0, formal=False)
            safe = dict(memory_available_bytes=20*1024**3, host_processes_observable=True,
                        gpu_observable=True, peers=[], processes=[])
            with patch('resources.run_guarded', side_effect=inline_guard), patch('resources.collect_snapshot', return_value=safe):
                self.invoke('benchmark', '--prepared', prepared, '--output', output, '--smoke-epochs', 1)
                self.assertEqual(len(calls), 60)
                self.invoke('benchmark', '--prepared', prepared, '--output', output, '--smoke-epochs', 1)
                self.assertEqual(len(calls), 60)
            ledger_path = output/'benchmark.json'
            original = ledger_path.read_bytes()
            ledger = json.loads(original)
            self.assertEqual({entry['status'] for entry in ledger['runs']}, {'completed'})
            self.assertEqual(len({workflow.digest(entry['identity']) for entry in ledger['runs']}), 60)
            summary_path = output/'summary.json'
            summary = summary_path.read_bytes()
            self.assertFalse(json.loads(summary)['formal_benchmark_complete'])
            self.assertTrue(json.loads(summary)['requested_scope_complete'])
            self.assertFalse(json.loads(summary)['formal_scope_complete'])
            with self.assertRaisesRegex(ValueError, 'Incomplete'):
                benchmark.evaluate_benchmark(output, require_complete=True)
            summary_path.unlink()
            (output/'report.md').unlink()
            with patch.object(training, 'train', side_effect=AssertionError('evaluation must not fit')):
                self.invoke('evaluate', '--benchmark', output)
                self.invoke('verify', '--benchmark', output, '--prepared', prepared)
            self.assertEqual(summary, summary_path.read_bytes())
            # A sealed ledger remains valid after the containing release is moved.
            ledger = json.loads(original)
            ledger['prepared_path'] = '../prepared/prepared.json'
            for entry in ledger['runs']:
                entry['run_path'] = 'runs/' + Path(entry['run_path']).name
            ledger_path.write_bytes(workflow.encoded(ledger) + b'\n')
            relocated = root/'relocated'
            output.rename(relocated)
            with patch.object(training, 'train', side_effect=AssertionError('evaluation must not fit')):
                self.invoke('evaluate', '--benchmark', relocated)
                self.invoke('verify', '--benchmark', relocated, '--prepared', prepared)
            self.assertEqual(summary, (relocated/'summary.json').read_bytes())
            self.assertEqual(workflow.encoded(ledger) + b'\n', (relocated/'benchmark.json').read_bytes())
            relocated.rename(output)
            ledger_path.write_bytes(original)
            for mutate in [lambda v: v['runs'].pop(), lambda v: v['runs'].append(v['runs'][0]),
                           lambda v: v['runs'][0]['identity'].update(seed=999),
                           lambda v: v.update(input_fingerprint='corrupt')]:
                ledger = json.loads(original)
                mutate(ledger)
                ledger_path.write_text(json.dumps(ledger))
                with self.assertRaises(ValueError):
                    benchmark.evaluate_benchmark(output)
            ledger_path.write_bytes(original)
            with self.assertRaisesRegex(ValueError, 'configuration'):
                benchmark.benchmark(prepared, output, smoke_epochs=2)
            prediction = Path(json.loads(original)['runs'][0]['run_path'])/'predictions.json'
            prediction.write_text('[]')
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                benchmark.evaluate_benchmark(output)

    def test_resource_deferral_failure_and_explicit_retry_preserve_evidence(self):
        import resources
        safe = dict(memory_available_bytes=20*1024**3, host_processes_observable=True,
                    gpu_observable=True, peers=[], processes=[])
        prepared = dict(expected_runs=[dict(condition='baseline', phenotype='cell_viability',
                        setting='cold_drug', fold=0, seed=17)], input_fingerprint='input', split_fingerprint='split')
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root/'benchmark'
            ledger_path = output/'benchmark.json'
            def ledger():
                return json.loads(ledger_path.read_text())
            def fail(command, log, **kwargs):
                entry = ledger()['runs'][0]
                run = Path(entry['run_path'])
                run.mkdir(parents=True)
                (run/'partial').write_text('failed checkpoint evidence')
                (log/'time.txt').write_text('failed measurement evidence')
                return dict(exit_code=1, stop_reason='monitor_abort', peak_tree_rss_bytes=123)
            with patch.object(training, '_prepared', return_value=prepared), \
                    patch.object(benchmark, 'evaluate_benchmark', return_value={'incomplete': True}), \
                    patch('resources.run_guarded', side_effect=fail) as launch:
                result = benchmark.benchmark(root/'prepared.json', output, smoke_epochs=1,
                    snapshot_fn=lambda **kw: dict(safe, memory_available_bytes=0))
                self.assertEqual(result, {'incomplete': True})
                launch.assert_not_called()
                entry = ledger()['runs'][0]
                self.assertEqual(entry['status'], 'pending')
                self.assertEqual(entry['blocked_reason'], 'available_ram_below_12_gib')
                self.assertFalse((output/'resources').exists())
                with self.assertRaisesRegex(ValueError, 'failed run'):
                    benchmark.benchmark(root/'prepared.json', output, smoke_epochs=1, snapshot_fn=lambda **kw: safe)
                entry = ledger()['runs'][0]
                self.assertEqual(entry['status'], 'failed')
                self.assertEqual(entry['resource_summary']['stop_reason'], 'monitor_abort')
                with self.assertRaisesRegex(ValueError, 'explicit retry'):
                    benchmark.benchmark(root/'prepared.json', output, smoke_epochs=1, snapshot_fn=lambda **kw: safe)
                self.assertEqual(launch.call_count, 1)
                with self.assertRaisesRegex(ValueError, 'failed run'):
                    benchmark.benchmark(root/'prepared.json', output, smoke_epochs=1,
                        snapshot_fn=lambda **kw: safe, retry_failed=True)
                archived = ledger()['runs'][0]['failed_attempts'][0]
                archive = Path(archived['archive_path'])
                self.assertEqual((archive/'run/partial').read_text(), 'failed checkpoint evidence')
                self.assertEqual((archive/'resources/time.txt').read_text(), 'failed measurement evidence')
                self.assertEqual(json.loads((archive/'attempt.json').read_text())['status'], 'failed')
                self.assertEqual(launch.call_count, 2)
                # A second guard preflight can refuse a launch after the first passes.
                launch.side_effect = resources.ResourceError('external_long_training')
                benchmark.benchmark(root/'prepared.json', output, smoke_epochs=1,
                    snapshot_fn=lambda **kw: safe, retry_failed=True)
                entry = ledger()['runs'][0]
                self.assertEqual(entry['status'], 'pending')
                self.assertEqual(entry['blocked_reason'], 'external_long_training')
                self.assertEqual(len(entry['failed_attempts']), 2)
            # Reporting preserves partial failed evidence and labels interrupted attempts explicitly.
            state = ledger()
            state['runs'][0]['status'] = 'running'
            run = Path(state['runs'][0]['run_path'])
            run.mkdir(parents=True, exist_ok=True)
            (run/'partial').write_text('interrupted evidence')
            ledger_path.write_text(json.dumps(state))
            with patch.object(training, '_prepared', return_value=prepared), \
                    patch('reporting.summarize', return_value={}) as summarize, \
                    patch('reporting.render_report', return_value='incomplete'):
                summary = benchmark.evaluate_benchmark(output)
                reported = summarize.call_args.args[1][0]
                self.assertEqual(reported['status'], 'failed')
                self.assertEqual(reported['error'], 'interrupted_run')
                self.assertFalse(summary['formal_benchmark_complete'])
                self.assertTrue((run/'partial').exists())
                self.assertEqual(ledger()['runs'][0]['status'], 'failed')

    def test_guarded_child_has_real_cap_threads_measurement_and_preflight_failure(self):
        # One genuine child validates the wrapper without loading any modeling dataset.
        import resources
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary)
            safe = dict(memory_available_bytes=20*1024**3, host_processes_observable=True,
                        gpu_observable=True, peers=[], processes=[])
            command = benchmark.capped_command([sys.executable, '-c',
                'import os,resource; assert resource.getrlimit(resource.RLIMIT_AS)[0] == 10*1024**3; '
                'assert os.environ["CUDA_VISIBLE_DEVICES"] == ""; assert os.environ["OMP_NUM_THREADS"] == "2"'], log)
            result = resources.run_guarded(command, log,
                                          snapshot_fn=lambda **kw: safe, interval=.01)
            self.assertEqual(result['exit_code'], 0)
            self.assertIn('Maximum resident set size', (log/'time.txt').read_text())
            with self.assertRaises(resources.ResourceError):
                resources.run_guarded([sys.executable, '-c', 'pass'], log/'blocked',
                                     snapshot_fn=lambda **kw: dict(safe, gpu_observable=False))


if __name__ == '__main__':
    unittest.main()
