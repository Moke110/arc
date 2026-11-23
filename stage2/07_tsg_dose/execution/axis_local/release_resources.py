"""Verify canonical RAM/GPU telemetry and the runner's service cap evidence."""
import csv
from datetime import datetime
import math
import re

from render_axis_local import METHODS, SOURCES

RESOURCE_FILES = (
    'formal-cap-start.txt', 'driver.time', 'resources-10s.log',
    'gpu-launch-window.log', 'gpu-10s.log',
)
OPTIONAL_RESOURCE_FILES = ('formal-cap-final.txt', 'formal-unit-journal.txt')
PHASES = (('validation', 'train'), ('validation', 'replay'),
          ('test', 'evaluate'), ('test', 'replay'))
GIB = 1024 ** 3


def samples(path, *, gpu):
    result = {}
    fields = {}
    stamp = None

    def save():
        if stamp is not None and {'MemTotal', 'MemAvailable', 'SwapFree'} <= fields.keys():
            if not gpu or 'gpu' in fields:
                result[stamp] = dict(fields)

    for line in path.read_text().splitlines():
        if re.fullmatch(r'\d{4}-\d{2}-\d{2}T\S+', line):
            save()
            stamp = datetime.fromisoformat(line).timestamp()
            fields = {}
            continue
        memory = re.fullmatch(r'(MemTotal|MemAvailable|SwapFree):\s+(\d+) kB', line)
        if memory:
            fields[memory[1]] = int(memory[2])
        gpu_values = re.fullmatch(r'(\d+),\s*(\d+),\s*(\d+)', line)
        if gpu_values:
            used, free, utilization = map(int, gpu_values.groups())
            if used + free > 0 and utilization <= 100:
                fields['gpu'] = (used, free, utilization)
    save()
    return result


def verify_window(series, begin, end, label):
    # systemd sampling sleeps ten seconds plus command execution time.
    stamps = sorted(stamp for stamp in series if begin - 15 <= stamp <= end + 15)
    if (not stamps or stamps[0] > begin + 15 or stamps[-1] < end - 15
            or any(right - left > 15 for left, right in zip(stamps, stamps[1:]))):
        raise ValueError(f'10-second {label} telemetry does not cover the formal phase window')
    return [series[stamp] for stamp in stamps]


def verify_resource_evidence(matrix, scenarios):
    resources = [matrix / name for name in RESOURCE_FILES]
    if any(path.is_symlink() or not path.is_file() for path in resources):
        raise ValueError('formal resource evidence is missing')
    cap, timing, ram, launch, gpu = resources
    settings = dict(line.split('=', 1) for line in cap.read_text().splitlines() if '=' in line)
    if (settings.get('Id') != 'arc-stage2-03-round1-formal.service'
            or not re.fullmatch(r'[0-9a-f]{32}', settings.get('InvocationID', ''))
            or settings.get('MemoryMax') != str(10 * GIB)
            or settings.get('MemorySwapMax') != str(8 * GIB)
            or not settings.get('ExecMainStartTimestamp')):
        raise ValueError('formal systemd cap or invocation identity missing')
    timing_text = timing.read_text()
    if (not re.search(r'^\s*Exit status:\s*0\s*$', timing_text, re.MULTILINE)
            or not re.search(r'^\s*Maximum resident set size \(kbytes\):\s*[1-9][0-9]*\s*$', timing_text, re.MULTILINE)):
        raise ValueError('successful formal driver wall/resource record required')
    bounds = []
    observed = set()
    for scenario in scenarios:
        path = matrix / 'reports' / scenario / 'fit_resources.csv'
        with path.open(newline='') as stream:
            for row in csv.DictReader(stream):
                if int(row['returncode']) != 0:
                    continue
                phase = (row['source'], row['scenario'], int(row['fold']), row['method'],
                         int(row['seed']), row['phase'], row['action'])
                if phase in observed or row['scenario'] != scenario:
                    raise ValueError('duplicate or cross-scenario formal phase resource record')
                observed.add(phase)
                start, elapsed = float(row['started_utc_unix_seconds']), float(row['seconds'])
                if not math.isfinite(start) or not math.isfinite(elapsed) or elapsed <= 0:
                    raise ValueError('invalid formal phase wall-time bounds')
                bounds.append((start, start + elapsed))
    expected = {(source, scenario, fold, method, 42, phase, action)
                for source in SOURCES for scenario in scenarios for fold in range(3)
                for method in METHODS for phase, action in PHASES}
    if observed != expected:
        raise ValueError('all 336 successful formal phase resource records required')
    begin, end = min(start for start, _ in bounds), max(finish for _, finish in bounds)
    journal = matrix / 'formal-unit-journal.txt'
    if journal.exists() and (journal.is_symlink() or not journal.stat().st_size):
        raise ValueError('invalid formal unit journal')
    ram_points = verify_window(samples(ram, gpu=False), begin, end, 'WSL RAM')
    gpu_samples = samples(launch, gpu=True) | samples(gpu, gpu=True)
    gpu_points = verify_window(gpu_samples, begin, end, 'GPU')
    return dict(cap_service=settings['Id'], invocation_id=settings['InvocationID'],
                memory_max_bytes=int(settings['MemoryMax']),
                memory_swap_max_bytes=int(settings['MemorySwapMax']),
                first_phase_unix_seconds=begin, last_phase_unix_seconds=end,
                ram_samples=len(ram_points), gpu_samples=len(gpu_points),
                min_available_ram_kib=min(point['MemAvailable'] for point in ram_points),
                max_gpu_used_mib=max(point['gpu'][0] for point in gpu_points))
