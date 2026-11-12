"""Write the auditable bounded-run acceptance record for #163.

The gate fails closed: the record is only written when the run artifacts agree
with the declaration (completed terminal status, declared finite prediction
shape over the full gene axis, case identities, zero reload error, the ADR 0063
early-stopping record, and measured resources).  Nothing is inferred from the
training log.
"""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from artifacts import read_json, sha256, write_json  # noqa: E402

TIME_PATTERN = {
    'elapsed': re.compile(r'Elapsed \(wall clock\) time.*?\):\s*(\S+)'),
    'maximum_resident_kbytes': re.compile(r'Maximum resident set size \(kbytes\):\s*(\d+)'),
    'exit_status': re.compile(r'Exit status:\s*(\d+)'),
}


def git(*args):
    return subprocess.run(['git', *args], capture_output=True, text=True,
                          check=True).stdout.strip()


def parse_time(path):
    text = Path(path).read_text()
    values = {}
    for name, pattern in TIME_PATTERN.items():
        match = pattern.search(text)
        values[name] = match.group(1) if match else None
    return values


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step', default=str(SCRIPT_DIR.parent.parent))
    parser.add_argument('--run', required=True)
    parser.add_argument('--release', required=True)
    parser.add_argument('--request', required=True)
    parser.add_argument('--config', required=True)
    parser.add_argument('--fit-time', required=True)
    parser.add_argument('--prepare-time', required=True)
    parser.add_argument('--prepare-output', required=True)
    parser.add_argument('--env-lock', required=True)
    parser.add_argument('--tests', help='JSON file holding the verification rows')
    parser.add_argument('--output', required=True)
    args = parser.parse_args(argv)

    run = Path(args.run)
    release = Path(args.release)
    step = Path(args.step)
    checkpoint = read_json(run / 'checkpoint.json')
    identity = read_json(run / 'identity.json')
    terminal = read_json(run / 'terminal.json')
    progress = read_json(run / 'progress.json')
    manifest = read_json(release / 'manifest.json')
    effective = read_json(run / 'config.json')
    request = read_json(args.request)
    tests = read_json(args.tests) if args.tests else []

    def require(condition, message):
        if not condition:
            raise SystemExit(f'bounded-run record gate failed: {message}')

    prediction = np.load(run / 'prediction.npy', mmap_mode='r')
    case_ids = read_json(run / 'case_ids.json')
    n_genes = int(manifest['n_genes'])
    require(terminal['status'] == 'complete', 'terminal status is not complete')
    require(progress['phase'] == 'complete', 'progress phase is not complete')
    require(prediction.shape == (len(case_ids), n_genes),
            f'prediction shape {prediction.shape} is not declared')
    require(prediction.dtype == np.float32, 'prediction is not float32')
    require(bool(np.isfinite(prediction).all()), 'prediction has non-finite values')
    require(len(case_ids) == manifest['roles']['test'], 'case count differs from the release')
    require(identity['prediction_sha256'] == sha256(run / 'prediction.npy'),
            'prediction digest differs from identity.json')
    require(identity['checkpoint_sha256'] == sha256(run / 'best.pt'),
            'checkpoint digest differs from identity.json')
    require(checkpoint['behavior']['reload_max_abs_error'] == 0.,
            'reload error is not exactly zero')
    require(checkpoint['early_stopping']['rules']['no_strict_decrease_epochs'] == 5,
            'the 5-epoch early-stopping rule is not recorded')
    require(checkpoint['early_stopping']['rules']['above_running_min_epochs'] == 20,
            'the 20-epoch early-stopping rule is not recorded')
    require('5-epoch rule dominates' in checkpoint['early_stopping']['dominance'],
            'the early-stopping dominance statement is missing')
    require(int(effective['max_epochs']) == int(checkpoint['config']['max_epochs']),
            'effective config and checkpoint config disagree')

    fit_time = parse_time(args.fit_time)
    prepare_time = parse_time(args.prepare_time)
    require(fit_time['maximum_resident_kbytes'] is not None,
            'the fit /usr/bin/time -v record has no peak RSS')
    environment = checkpoint['environment']
    files = {
        'release/manifest.json': sha256(release / 'manifest.json'),
        'run/prediction.npy': sha256(run / 'prediction.npy'),
        'run/case_ids.json': sha256(run / 'case_ids.json'),
        'run/best.pt': sha256(run / 'best.pt'),
        'run/checkpoint.json': sha256(run / 'checkpoint.json'),
        'run/identity.json': sha256(run / 'identity.json'),
        'run/terminal.json': sha256(run / 'terminal.json'),
        'run/config.json': sha256(run / 'config.json'),
        'request.json': sha256(args.request),
        'config/bounded_run.json': sha256(args.config),
    }
    record = {
        'schema': 'stage2-02-models-tenv-bounded-run-v1',
        'ticket': 163,
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'branch': git('rev-parse', '--abbrev-ref', 'HEAD'),
        'head_sha': git('rev-parse', 'HEAD'),
        'status': terminal['status'],
        'identity': checkpoint['identity'],
        'declared': {
            'bound_epochs': effective['max_epochs'],
            'full_epoch_ceiling': checkpoint['early_stopping']['rules']['max_epochs'],
            'batch_size': effective['batch_size'],
            'wells': manifest['roles'],
            'n_genes': n_genes,
            'early_stopping': checkpoint['early_stopping']['rules'],
            'early_stopping_dominance': checkpoint['early_stopping']['dominance'],
        },
        'outcome': {
            'best_epoch': checkpoint['best_epoch'],
            'best_validation_loss': checkpoint['best_validation_loss'],
            'stop_reason': checkpoint['stop_reason'],
            'history_epochs': len(checkpoint['history']),
            'training_seconds': checkpoint['training_seconds'],
            'inference_seconds': checkpoint['inference_seconds'],
            'behavior': checkpoint['behavior'],
            'selection': checkpoint['selection'],
            'output': checkpoint['output'],
        },
        'prediction': {
            'shape': list(prediction.shape),
            'dtype': str(prediction.dtype),
            'finite': True,
            'case_ids': len(case_ids),
            'first_case_id': case_ids[0],
            'last_case_id': case_ids[-1],
            'sha256': files['run/prediction.npy'],
        },
        'resources': {
            'fit_time': fit_time,
            'prepare_time': prepare_time,
            'checkpoint_peak_rss_kbytes': checkpoint['peak_rss_kbytes'],
            'checkpoint_peak_gpu_bytes': checkpoint['peak_gpu_bytes'],
            'threads': effective['threads'],
            'prepare_peak_rss_kbytes': prepare_time['maximum_resident_kbytes'],
        },
        'environment': {
            'lock': str(Path(args.env_lock)),
            'lock_sha256': sha256(args.env_lock),
            'interpreter': environment['executable'],
            'python': environment['python'],
            'torch': environment['torch'],
            'numpy': environment['numpy'],
            'device': environment['device'],
            'cuda_available': environment['cuda_available'],
            'threads': environment['threads'],
            'numba_cache_dir': environment['numba_cache_dir'],
            'upstream_source': environment['upstream_source'],
            'upstream_prnet_module_sha256': environment['upstream_prnet_module_sha256'],
            'prepare_interpreter': read_json(step / 'log/env.prepare.json')['executable'],
        },
        'deviations': [
            {'deviation': 'no GPU this boot; CPU fallback',
             'evidence': {'torch_cuda_available': environment['cuda_available'],
                          'dev_dxg_present': Path('/dev/dxg').exists(),
                          'device': environment['device']}},
            {'deviation': 'NUMBA_CACHE_DIR required for scanpy/numba import in the '
                          'pinned environment',
             'evidence': environment['numba_cache_dir']},
            {'deviation': 'bounded run: 3 epochs on a declared well prefix, not the '
                          '2f x 2s matrix (#165)',
             'evidence': {'epochs': effective['max_epochs'],
                          'wells': manifest['roles']}},
            {'deviation': 'worktree lacks ignored external source trees; the fit '
                          'command resolved the primary checkout tree',
             'evidence': environment['upstream_source']},
        ],
        'paths': {
            'run': str(run),
            'release': str(release),
            'request': str(args.request),
            'bounded_config': str(args.config),
            'prepare_output': str(Path(args.prepare_output)),
        },
        'files_sha256': files,
        'release_source_sha256': manifest.get('source'),
        'verification': tests,
    }
    write_json(args.output, record)
    print(json.dumps({key: record[key] for key in
                      ('status', 'identity', 'declared', 'resources')},
                     sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
