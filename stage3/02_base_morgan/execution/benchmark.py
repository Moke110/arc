"""Sequential, resumable benchmark accounting with immutable verified runs."""
import json
import os
from pathlib import Path
import platform
import re
import sys
from workflow import digest, encoded


def _write(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_bytes(encoded(value) + b'\n')
    temporary.replace(path)


def _identity_key(identity):
    return digest(identity)


def _validate_ledger(prepared, ledger, output):
    if ledger.get('schema_version') != 1 or ledger.get('mode') not in ('formal', 'smoke'):
        raise ValueError('Invalid benchmark execution configuration')
    epochs = ledger.get('smoke_epochs')
    if ((ledger['mode'] == 'formal' and epochs is not None)
            or (ledger['mode'] == 'smoke' and (type(epochs) is not int or not 1 <= epochs <= 300))):
        raise ValueError('Invalid benchmark execution configuration')
    expected = prepared['expected_runs']
    entries = ledger['runs']
    keys = [_identity_key(entry['identity']) for entry in entries]
    if len(keys) != len(set(keys)):
        raise ValueError('Duplicate benchmark run identity')
    if set(keys) != {_identity_key(identity) for identity in expected}:
        raise ValueError('Missing or unexpected benchmark run identities')
    if (ledger['input_fingerprint'] != prepared['input_fingerprint']
            or ledger['split_fingerprint'] != prepared['split_fingerprint']
            or ledger['prepared_fingerprint'] != digest(prepared)):
        raise ValueError('Incompatible benchmark prepared fingerprint')
    for entry in entries:
        if entry['status'] not in ['pending', 'running', 'completed', 'failed']:
            raise ValueError('Invalid run status')
        if entry['run_path'] != str((output/'runs'/_identity_key(entry['identity'])).resolve()):
            raise ValueError('Unexpected benchmark run path')
    actual = {str(path.resolve()) for path in (output/'runs').glob('*')}
    allowed = {entry['run_path'] for entry in entries}
    if actual - allowed:
        raise ValueError('Unexpected run directories')


def _completed(entry, prepared_path, ledger):
    from training import verify
    path = Path(entry['run_path'])
    verify(path, prepared_path)
    manifest = json.loads((path/'run.json').read_text())
    if (manifest['identity'] != entry['identity'] or manifest['mode'] != ledger['mode']
            or manifest['smoke_epochs'] != ledger['smoke_epochs']):
        raise ValueError('Incompatible existing run configuration')


def evaluate_benchmark(output, prepared_path=None, require_complete=False):
    """Rebuild report from verified predictions; never launch training."""
    from training import _prepared
    from reporting import summarize, render_report
    output = Path(output).resolve()
    ledger = json.loads((output/'benchmark.json').read_text())
    prepared_path = Path(prepared_path or ledger['prepared_path'])
    prepared = _prepared(prepared_path)
    _validate_ledger(prepared, ledger, output)
    for entry in ledger['runs']:
        if entry['status'] == 'completed':
            _completed(entry, prepared_path, ledger)
        elif entry['status'] == 'pending' and Path(entry['run_path']).exists():
            raise ValueError('Pending run has unexpected artifacts')
        elif entry['status'] == 'running':
            entry['status'], entry['error'] = 'failed', 'interrupted_run'
    _write(output/'benchmark.json', ledger)
    summary = summarize(prepared, ledger['runs'], require_complete=require_complete)
    summary['benchmark_mode'] = ledger['mode']
    summary['formal_benchmark_complete'] = ledger['mode'] == 'formal' and summary.get('complete', False)
    _write(output/'summary.json', summary)
    (output/'report.md').write_text(render_report(summary))
    return summary


def _archive_attempt(entry, output, log):
    """Move failed evidence aside only on an explicitly requested retry."""
    attempts = entry.setdefault('failed_attempts', [])
    index = len(attempts) + 1
    archive = output/'failed_attempts'/_identity_key(entry['identity'])/str(index)
    if archive.exists():
        raise ValueError('Failed attempt archive already exists')
    archive.mkdir(parents=True)
    metadata = {key: value for key, value in entry.items() if key != 'failed_attempts'}
    metadata['attempt_index'] = index
    _write(archive/'attempt.json', metadata)
    for path, name in [(Path(entry['run_path']), 'run'), (log, 'resources')]:
        if path.exists():
            path.rename(archive/name)
    attempts.append(dict(archive_path=str(archive), **metadata))
    for key in ('error', 'resource_summary', 'peak_rss_bytes', 'blocked_reason'):
        entry.pop(key, None)
    entry['status'] = 'pending'


def capped_command(command, log):
    # Apply portable hard address-space cap before imports; env is child-local.
    wrapper = ('import os,resource,sys; '
               'resource.setrlimit(resource.RLIMIT_AS,(10*1024**3,10*1024**3)); '
               'os.environ.update(OMP_NUM_THREADS="2",MKL_NUM_THREADS="2",'
               'OPENBLAS_NUM_THREADS="2",NUMEXPR_NUM_THREADS="2",CUDA_VISIBLE_DEVICES=""); '
               'os.execv(sys.argv[1],sys.argv[1:])')
    return [sys.executable, '-c', wrapper, '/usr/bin/time', '-v', '-o', str(log/'time.txt'), *command]


def benchmark(prepared_path, output, smoke_epochs=None, snapshot_fn=None, interval=5, retry_failed=False):
    """Execute one capped child per run. snapshot_fn is a test-only injection seam."""
    from training import _prepared
    from resources import run_guarded, ResourceError, collect_snapshot, preflight
    prepared_path, output = Path(prepared_path).resolve(), Path(output).resolve()
    prepared = _prepared(prepared_path)
    if smoke_epochs is not None and (type(smoke_epochs) is not int or not 1 <= smoke_epochs <= 300):
        raise ValueError('smoke_epochs must be between 1 and 300')
    mode = 'formal' if smoke_epochs is None else 'smoke'
    ledger_path = output/'benchmark.json'
    if output.exists():
        if not ledger_path.exists():
            raise ValueError('Output exists with incompatible benchmark artifacts')
        ledger = json.loads(ledger_path.read_text())
        _validate_ledger(prepared, ledger, output)
        if ledger['mode'] != mode or ledger['smoke_epochs'] != smoke_epochs:
            raise ValueError('Incompatible benchmark execution configuration')
    else:
        output.mkdir(parents=True)
        ledger = dict(schema_version=1, mode=mode, smoke_epochs=smoke_epochs,
                      prepared_path=str(prepared_path), prepared_fingerprint=digest(prepared),
                      input_fingerprint=prepared['input_fingerprint'], split_fingerprint=prepared['split_fingerprint'],
                      execution=dict(python=platform.python_version(), platform=platform.platform(),
                          cpu_count=os.cpu_count(), device='cpu', threads=2, cuda_visible_devices='',
                          memory_cap_bytes=10*1024**3, cap_mechanism='RLIMIT_AS', sequential=True),
                      runs=[dict(identity=identity, status='pending',
                                 run_path=str(output/'runs'/_identity_key(identity)))
                            for identity in prepared['expected_runs']])
        _write(ledger_path, ledger)
    # Validate ALL existing artifacts before launching another child.
    for entry in ledger['runs']:
        if entry['status'] == 'completed':
            _completed(entry, prepared_path, ledger)
        elif entry['status'] == 'pending' and Path(entry['run_path']).exists():
            raise ValueError('Pending run has unexpected artifacts')
        elif entry['status'] == 'running':
            entry['status'], entry['error'] = 'failed', 'interrupted_run'
    _write(ledger_path, ledger)
    if not retry_failed and any(entry['status'] == 'failed' for entry in ledger['runs']):
        raise ValueError('Failed runs require explicit retry_failed=True after investigation')
    cli = Path(__file__).with_name('workflow.py')
    for entry in ledger['runs']:
        if entry['status'] == 'completed':
            continue
        identity = entry['identity']
        try:
            preflight((snapshot_fn or collect_snapshot)(exclude_pids=()), formal=mode == 'formal')
        except ResourceError as error:
            entry['blocked_reason'] = str(error)
            _write(ledger_path, ledger)
            return evaluate_benchmark(output, prepared_path)
        log = output/'resources'/_identity_key(identity)
        if entry['status'] == 'failed':
            _archive_attempt(entry, output, log)
        log.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, str(cli), 'train', '--prepared', str(prepared_path),
                   '--output', entry['run_path']]
        for key, value in identity.items():
            command.extend(['--'+key, str(value)])
        if smoke_epochs is not None:
            command.extend(['--smoke-epochs', str(smoke_epochs)])
        capped = capped_command(command, log)
        entry['status'] = 'running'
        entry.pop('error', None)
        entry.pop('blocked_reason', None)
        _write(ledger_path, ledger)
        try:
            result = run_guarded(capped, log, formal=mode == 'formal',
                                 snapshot_fn=snapshot_fn, interval=interval)
            entry['resource_summary'] = result
            if result['exit_code'] != 0 or result['stop_reason']:
                raise ValueError(f'Guarded child failed: {result}')
            _completed(entry, prepared_path, ledger)
            peak = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', (log/'time.txt').read_text())
            if peak is None:
                raise ValueError('Missing measured peak RSS')
            entry['peak_rss_bytes'] = int(peak[1])*1024
            entry['status'] = 'completed'
        except ResourceError as error:
            # run_guarded can refuse its second preflight before starting a child.
            entry['status'], entry['blocked_reason'] = 'pending', str(error)
            _write(ledger_path, ledger)
            return evaluate_benchmark(output, prepared_path)
        except (ValueError, OSError) as error:
            entry['status'], entry['error'] = 'failed', str(error)
            _write(ledger_path, ledger)
            raise ValueError(f'Benchmark stopped; failed run {identity}: {error}') from error
        _write(ledger_path, ledger)
    return evaluate_benchmark(output, prepared_path)
