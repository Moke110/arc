"""Read the formal run's recorded caps, timing and ten-second telemetry."""
from datetime import datetime
import hashlib
import json
import re


def key_values(path):
    return dict(line.split('=', 1) for line in path.read_text().splitlines() if '=' in line)


def timed_value(timing, label):
    match = re.search(rf'^\s*{re.escape(label)}:\s*(\S.*)$', timing, re.MULTILINE)
    if match is None:
        raise ValueError(f'GNU time result missing {label}')
    return match.group(1)


def telemetry(path, start, end):
    sampled = {}
    current = None
    for line in path.read_text().splitlines():
        if re.match(r'^\d{4}-\d\d-\d\dT', line):
            current = datetime.fromisoformat(line).timestamp()
            sampled[current] = {}
        elif current is not None:
            ram = re.fullmatch(r'MemAvailable:\s+(\d+) kB', line)
            gpu = re.fullmatch(r'(\d+),\s*(\d+),\s*(\d+)', line)
            if ram:
                sampled[current]['available_ram_kib'] = int(ram[1])
            elif gpu:
                used, free, utilization = (int(value) for value in gpu.groups())
                if used + free > 0 and 0 <= utilization <= 100:
                    sampled[current]['gpu_used_mib'] = used
    points = [point for point, values in sampled.items() if len(values) == 2]
    if (not points or points[0] > start + 15 or points[-1] < end - 15
            or any(not 0 < right - left <= 15 for left, right in zip(points, points[1:]))):
        raise ValueError(f'fit lacks continuous ten-second WSL/GPU telemetry: {path}')
    return dict(samples=len(points), min_ram_kib=min(sampled[point]['available_ram_kib'] for point in points),
                max_gpu_mib=max(sampled[point]['gpu_used_mib'] for point in points))


def verify_resources(root, evidence):
    fits = {(row['kind'], row['fold']): row for row in evidence['fits']}
    phases = {
        'tahoe': ('fit', 'replay'),
        'tenv': ('fit', 'replay', 'evaluate', 'replay-test'),
        'transfer': ('fit', 'replay', 'evaluate', 'replay-test'),
        'tahoe-only': ('evaluate', 'replay-test'),
    }
    records = []
    for step, actions in phases.items():
        for fold in range(3):
            for action in actions:
                label = f'{step}-formal-{action}-{fold}'
                attempts = [
                    path for path in root.glob(label + '-*')
                    if (path / 'terminal.txt').exists()
                    and key_values(path / 'terminal.txt').get('status') == '0'
                ]
                if len(attempts) != 1:
                    raise ValueError(f'require exactly one successful run for {label}')
                attempt = attempts[0]
                preflight = key_values(attempt / 'preflight.txt')
                terminal = key_values(attempt / 'terminal.txt')
                cap = key_values(attempt / 'scope.txt')
                start = datetime.fromisoformat(preflight['started']).timestamp()
                end = datetime.fromisoformat(terminal['ended']).timestamp()
                timing = (attempt / 'time.txt').read_text()
                peak = int(timed_value(timing, 'Maximum resident set size (kbytes)'))
                elapsed = timed_value(timing, 'Elapsed (wall clock) time (h:mm:ss or m:ss)')
                if (end < start or int(preflight['ram_kib']) < 12 * 1024 ** 2
                        or int(preflight['gpu_mib']) <= 0 or peak <= 0
                        or timed_value(timing, 'Exit status') != '0'
                        or cap['Id'] != f'arc-stage2-transfer-{attempt.name}.scope'
                        or not re.fullmatch('[0-9a-f]{32}', cap['InvocationID'])
                        or int(cap['MemoryMax']) != 10 * 1024 ** 3
                        or int(cap['MemorySwapMax']) != 8 * 1024 ** 3):
                    raise ValueError(f'resource admission, scope or timing is invalid for {label}')
                record = dict(phase=label, start=preflight['started'], end=terminal['ended'],
                              peak_rss_kib=peak, wall_clock=elapsed, cap=cap,
                              preflight_ram_kib=int(preflight['ram_kib']),
                              preflight_gpu_mib=int(preflight['gpu_mib']))
                files = ['preflight.txt', 'terminal.txt', 'scope.txt', 'time.txt']
                if action == 'fit':
                    info = json.loads((attempt / 'stdout.txt').read_text())
                    fit = fits['tenv-only' if step == 'tenv' else step, fold]
                    if any(info[key] != fit[key] for key in (
                        'epochs', 'best_epoch', 'training_steps', 'validation_mse',
                        'elapsed_seconds', 'peak_rss_kib',
                    )):
                        raise ValueError(f'formal timing and matrix fit differ for {label}')
                    samples = attempt / 'resources-10s.log'
                    record['telemetry'] = telemetry(samples, start, end)
                    files.extend(('resources-10s.log', 'stdout.txt'))
                record['evidence_sha256'] = {
                    name: hashlib.sha256((attempt / name).read_bytes()).hexdigest()
                    for name in files
                }
                records.append((start, end, record))
    records.sort(key=lambda item: (item[0], item[1], item[2]['phase']))
    if any(previous[1] > following[0] for previous, following in zip(records, records[1:])):
        raise ValueError('formal resource windows overlap')
    return [record for _, _, record in records]
