"""Read the formal run's recorded caps, timing and ten-second telemetry."""
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import shlex


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



def verify_phase_output(attempt, candidate, step, action, fold):
    """Bind the executed command and printed result to this formal artifact."""
    step_script = 'transfer' if step == 'tahoe-only' else step
    command = shlex.split((attempt / 'command.txt').read_text())
    if len(command) < 3 or Path(command[1]).parts[-2:] != ('execution', f'{step_script}.py'):
        raise ValueError('formal phase command uses the wrong entry point')
    if command[2] != action or len(command[3:]) % 2:
        raise ValueError('formal phase command has the wrong action or flag layout')
    flags = dict(zip(command[3::2], command[4::2], strict=True))
    if (len(flags) * 2 != len(command[3:]) or flags.get('--fold') != str(fold)
            or flags.get('--device', 'cuda') != 'cuda'
            or Path(flags.get('--candidate', '')).name == ''
            or any(flag.startswith('--smoke') for flag in flags)):
        raise ValueError('formal phase command has the wrong fold, candidate, or device')
    if flags.get('--path', 'transfer') != step and (
        step in ('transfer', 'tahoe-only') or '--path' in flags
    ):
        raise ValueError('formal phase command has the wrong path')
    payload = json.loads((attempt / 'stdout.txt').read_text())
    model_path = 'tenv-only' if step == 'tenv' else step
    relative = Path(model_path) / f'fold-{fold}'
    if action == 'fit':
        return payload
    if action == 'replay':
        if step != 'tahoe':
            relative /= 'fit'
        expected = json.loads((candidate / relative / 'validation-replay.json').read_text())
    else:
        manifest = relative / 'test/manifest.json'
        # Evaluation prints its artifact path. The separate replay phase prints
        # the sealed manifest, and matrix verification checks its replay proof.
        if action == 'evaluate':
            exported = Path(flags['--candidate']) / manifest
            if payload != str(exported):
                raise ValueError('test export does not name its expected manifest')
            return payload
        expected = json.loads((candidate / manifest).read_text())
    if payload != expected:
        raise ValueError('formal phase stdout differs from its verified receipt')
    return payload


def verify_resources(root, evidence, candidate):
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
                files = ['preflight.txt', 'terminal.txt', 'scope.txt', 'time.txt', 'command.txt', 'stdout.txt']
                info = verify_phase_output(attempt, candidate, step, action, fold)
                if action == 'fit':
                    fit = fits['tenv-only' if step == 'tenv' else step, fold]
                    if any(info[key] != fit[key] for key in (
                        'epochs', 'best_epoch', 'training_steps', 'validation_mse',
                        'elapsed_seconds', 'peak_rss_kib',
                    )):
                        raise ValueError(f'formal timing and matrix fit differ for {label}')
                    samples = attempt / 'resources-10s.log'
                    record['telemetry'] = telemetry(samples, start, end)
                    files.append('resources-10s.log')
                record['evidence_sha256'] = {
                    name: hashlib.sha256((attempt / name).read_bytes()).hexdigest()
                    for name in files
                }
                records.append((start, end, record))
    records.sort(key=lambda item: (item[0], item[1], item[2]['phase']))
    if any(previous[1] > following[0] for previous, following in zip(records, records[1:])):
        raise ValueError('formal resource windows overlap')
    return [record for _, _, record in records]
