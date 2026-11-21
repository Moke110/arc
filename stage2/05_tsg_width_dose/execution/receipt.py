"""Assemble a compact Step 05 report acceptance receipt from verified artifacts."""

import argparse
from datetime import datetime
from hashlib import sha256
import json
from pathlib import Path
import re

from run import ACCEPTED, SCENARIO, SEED, digest, plan, run_path

FIGURES = ('learning_curves', 'metrics_balloon', 'metrics_bar', 'metrics_per_dose_ballone')
MEMORY_MAX = 10 * 1024**3
SWAP_MAX = 8 * 1024**3


def inventory(root):
    return {str(path.relative_to(root)): digest(path)
            for path in sorted(root.rglob('*')) if path.is_file()}


def kv(path):
    return dict(line.split('=', 1) for line in path.read_text().splitlines())


def phase_record(path):
    files = inventory(path)
    scope = kv(path / 'scope.txt')
    preflight = kv(path / 'preflight.txt')
    terminal = kv(path / 'terminal.txt')
    if (int(scope['MemoryMax']), int(scope['MemorySwapMax'])) != (MEMORY_MAX, SWAP_MAX):
        raise ValueError(f'Unexpected memory caps in {path}')
    if not re.fullmatch('[0-9a-f]{32}', scope['InvocationID']):
        raise ValueError(f'Invalid systemd invocation ID in {path}')
    start = datetime.fromisoformat(preflight['started'])
    end = datetime.fromisoformat(terminal['ended'])
    lines = (path / 'resources-10s.log').read_text().splitlines()
    if not lines or len(lines) % 5:
        raise ValueError(f'Invalid telemetry series in {path}')
    samples = []
    for date, total, ram, swap, gpu in zip(*(iter(lines),) * 5, strict=True):
        timestamp = datetime.fromisoformat(date)
        total_kib = int(total.split()[1])
        ram_kib = int(ram.split()[1])
        swap_kib = int(swap.split()[1])
        gpu_used, gpu_free, gpu_util = map(int, gpu.split(','))
        samples.append((timestamp, total_kib, ram_kib, swap_kib, gpu_used, gpu_free, gpu_util))
    intervals = [(next_[0] - prev[0]).total_seconds() for prev, next_ in zip(samples, samples[1:])]
    if not start <= samples[0][0] <= samples[-1][0] <= end or any(gap <= 0 for gap in intervals):
        raise ValueError(f'Telemetry outside attempt in {path}')
    match = re.search(r'Maximum resident set size \(kbytes\): (\d+)', (path / 'time.txt').read_text())
    status = int(terminal['status'])
    if status == 0 and (not match or (path / 'monitor-failure.txt').exists()):
        raise ValueError(f'Success lacks valid timing and monitoring evidence in {path}')
    attempted = [json.loads(line) for line in (path / 'stdout.txt').read_text().splitlines() if line]
    if not attempted:
        raise ValueError(f'No coordinate registry in {path}')
    return dict(attempted=attempted, status=status, started=preflight['started'], ended=terminal['ended'],
                memory_max_bytes=MEMORY_MAX, swap_max_bytes=SWAP_MAX, scope=scope['Id'],
                invocation_id=scope['InvocationID'], files=files,
                telemetry_samples=len(samples), max_sample_gap_seconds=max(intervals, default=None),
                wsl_total_ram_kib=samples[0][1], min_available_ram_kib=min(r[2] for r in samples),
                min_free_swap_kib=min(r[3] for r in samples), max_gpu_used_mib=max(r[4] for r in samples),
                peak_rss_kib=int(match.group(1)) if match else None)


def retain_review_annotations(receipt, path):
    reviewed = json.loads(path.read_text())
    for key in ('report_phase_peak_rss_kib', 'visual_inspection'):
        receipt[key] = reviewed[key]
    for name, phase in receipt['resources'].items():
        published = reviewed['resources'][name]['published_telemetry']
        if Path(published).name != published or (
            digest(path.parent / published) != phase['files']['resources-10s.log']
        ):
            raise ValueError(f'Published telemetry does not match capped attempt {name}')
        phase['published_telemetry'] = published
    if receipt != reviewed:
        raise ValueError(f'Reviewed receipt differs from the reconstructed evidence: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('candidate', 'verified-plan', 'science-evidence',
                'report', 'replay', 'script-dir', 'output'):
        parser.add_argument('--' + arg, type=Path, required=True)
    parser.add_argument('--reviewed-receipt', type=Path,
                        help='Check against a published receipt and retain measured and visual annotations')
    args = parser.parse_args()
    verified = json.loads(args.verified_plan.read_text())
    if verified != plan():
        raise ValueError('Verified coordinates differ from the registered plan')
    candidate_plan = args.candidate / 'plan.json'
    identity = json.loads(candidate_plan.read_text())
    input_digest = digest(ACCEPTED)
    if identity['input_manifest_sha256'] != input_digest:
        raise ValueError('Candidate differs from frozen input')
    scientific = json.loads(args.science_evidence.read_text())
    comparison = scientific['original_comparison']
    if len(comparison) != 12:
        raise ValueError('Expected twelve original-width Step 04 comparisons')
    comparison_bytes = json.dumps(comparison, sort_keys=True, allow_nan=False).encode()
    phases = {path.name: phase_record(path) for path in sorted((args.candidate / 'phases').iterdir())}
    if not phases or all(row['status'] for row in phases.values()):
        raise ValueError('No successful capped execution')
    origins = []
    for origin in verified:
        path = run_path(args.candidate, origin)
        terminal = json.loads((path / 'terminal.json').read_text())
        if terminal['status'] != 'complete' or terminal['validation_replay'] != 'exact' or (
            terminal['training_pid'] == terminal['evaluation_pid']
        ):
            raise ValueError(f'Invalid fresh-process replay at {path}')
        manifest = path / 'manifest.json'
        sealed_at = manifest.stat().st_mtime
        observed_attempts = [name for name, phase in phases.items() if origin in phase['attempted']
                             and datetime.fromisoformat(phase['started']).timestamp() <= sealed_at
                             < datetime.fromisoformat(phase['ended']).timestamp() + 1]
        if len(observed_attempts) != 1:
            raise ValueError(f'No unique capped execution window for {manifest}')
        origins.append(dict(origin=origin, attempt=observed_attempts[0], manifest_sha256=digest(manifest),
                            test_samples=terminal['test_samples'], best_epoch=terminal['best_epoch'],
                            epochs=terminal['epochs'], validation_mse=terminal['validation_mse'],
                            validation_replay='exact'))
    plots = {f'{name}.png': digest(args.report / f'{name}.png') for name in FIGURES}
    if plots != inventory(args.replay):
        raise ValueError('CSV-only replay inventory does not match publication figures')
    tables = inventory(args.report / 'data')
    if set(tables) != {f'{name}.csv' for name in FIGURES}:
        raise ValueError('Unexpected figure data inventory')
    receipt = dict(schema='tsg-width-dose-report-v1', status='passed', source='tenv',
                   scenario=SCENARIO, seed=SEED, formal_fits=len(verified),
                   input_manifest_sha256=input_digest, plan_sha256=digest(candidate_plan),
                   origins=origins, resources=phases, figures=plots, csv_only_replay=True,
                   plotted_tables=tables,
                   scientific_verification_sha256=digest(args.science_evidence),
                   original_width_comparison_sha256=sha256(comparison_bytes).hexdigest(),
                   report_scripts_sha256={name: digest(args.script_dir / name) for name in
                                          ('plot_contract.py', 'prepare.py', 'render.py')})
    if args.reviewed_receipt:
        retain_review_annotations(receipt, args.reviewed_receipt)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
