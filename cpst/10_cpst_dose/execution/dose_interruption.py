"""Verify observed resource evidence for the interrupted Step 10 matrix."""

from dataclasses import dataclass
from datetime import datetime
from hashlib import file_digest
import itertools
import json
import math
from pathlib import Path
import re

STARTED = '39f53479d3a045ac8e11786248231fbf'
ACCOUNTED = 'ae8f7b866b0347b9af31fe1c80b127c0'
SCHEMA = 'cpst-step10-matrix-interruption-v1'


@dataclass(frozen=True)
class Scope:
    session: str
    start: float
    stop: float
    unix_start: float
    unix_stop: float


def digest(path):
    with path.open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def json_lines(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def pinned_file(root, name, pin):
    path = root / name
    if Path(name).name != name or path.is_symlink() or digest(path) != pin:
        raise ValueError(f'interruption journal digest or path differs: {name}')
    return path


def scope_windows(root, journals):
    scopes, boots, commands, identities = [], set(), set(), set()
    if len(journals) != 2:
        raise ValueError('interrupted matrix requires exactly two scope journals')
    for entry in journals:
        unit = entry['unit']
        if not re.fullmatch(r'run-[a-zA-Z0-9]+\.scope', unit):
            raise ValueError('invalid matrix systemd scope')
        path = pinned_file(root, entry['path'], entry['sha256'])
        rows = json_lines(path)
        events = [row for row in rows if row.get('USER_UNIT') == unit
                  and row.get('MESSAGE_ID') in (STARTED, ACCOUNTED)]
        if [row['MESSAGE_ID'] for row in events] != [STARTED, ACCOUNTED]:
            raise ValueError(f'expected one unit start and terminal accounting record: {unit}')
        boots.update(row['_BOOT_ID'] for row in events)
        start, end = events
        identity = start['USER_INVOCATION_ID']
        if (identity != end['USER_INVOCATION_ID'] or identity in identities
                or int(end['CPU_USAGE_NSEC']) <= 0
                or not start['MESSAGE'].startswith(f'Started {unit} - ')):
            raise ValueError(f'scope journal invocation or terminal accounting differs: {unit}')
        identities.add(identity)
        command = start['MESSAGE'].removeprefix(f'Started {unit} - ')
        match = re.search(r'\bdose_matrix\.py\s+(.*)\.$', command)
        if match is None or '--output' not in command or '--logs' not in command:
            raise ValueError(f'scope journal does not identify the matrix worker: {unit}')
        commands.add(match[1])
        scopes.append(Scope(entry['monitor_session'], *(
            int(row[key]) / 1_000_000
            for key in ('__MONOTONIC_TIMESTAMP', '__REALTIME_TIMESTAMP') for row in events
        )))
    first, last = scopes
    if (len(boots) != 1 or len(commands) != 1 or first.session == last.session
            or not first.start < first.stop < last.start < last.stop
            or not first.unix_start < first.unix_stop < last.unix_start < last.unix_stop):
        raise ValueError('scope journals do not prove ordered, idle-separated same-boot attempts')
    return scopes


def active_samples(supervisor, matrix, scopes, final_measurement):
    pools = [list(group) for _, group in itertools.groupby(matrix, key=lambda row: row['monitor_session'])]
    if [pool[0]['monitor_session'] for pool in pools] != [scope.session for scope in scopes]:
        raise ValueError('matrix monitor sessions do not match the authorized scopes')
    if final_measurement < scopes[-1].unix_start or final_measurement > scopes[-1].unix_stop:
        raise ValueError('kernel peak receipt falls outside the completed scope')
    selected, evidence = [], []
    for index, (scope, pool) in enumerate(zip(scopes, pools, strict=True)):
        ticks = [row['monotonic_seconds'] for row in pool]
        if (not scope.start <= ticks[0] <= scope.start + 10
                or not scope.stop - 12 <= ticks[-1] <= scope.stop
                or any(not math.isfinite(right - left) or right <= left
                       for left, right in zip(ticks, ticks[1:]))):
            raise ValueError('matrix observations do not cover an evidenced active attempt')
        before = [row for row in supervisor if row['monotonic_seconds'] <= ticks[0]]
        active = [row for row in supervisor if ticks[0] < row['monotonic_seconds'] <= scope.stop]
        if not before or not active:
            raise ValueError('missing independent observations during an active attempt')
        rows = [before[-1], *active]
        if index == len(scopes) - 1:
            # Keep the existing successful-matrix bracket through its final measurement.
            while rows[-1]['monotonic_seconds'] < ticks[-1] or (
                datetime.fromisoformat(rows[-1]['wall_time']).timestamp() < final_measurement
            ):
                after = [row for row in supervisor if row['monotonic_seconds'] > rows[-1]['monotonic_seconds']]
                if not after:
                    raise ValueError('missing independent observations through successful matrix completion')
                rows.append(after[0])
        times = [row['monotonic_seconds'] for row in rows]
        if (not scope.start - 12 <= times[0] <= scope.start + 10
                or times[-1] < scope.stop - 12
                or any(not math.isfinite(right - left) or abs(right - left - 10) > 1
                       for left, right in zip(times, times[1:]))):
            raise ValueError('independent ten-second cadence is incomplete inside an active scope')
        for row in rows:
            total, used, _ = map(int, row['gpu_total_used_mib_and_utilization_percent'].split(','))
            if (row['wsl']['MemAvailable_kib'] < 12 * 1024**2 or total - used < 512
                    or int(row['windows_host_total_and_available_kib'].split(',')[1]) <= 0):
                raise ValueError('independent resource sample exceeds the formal envelope')
        evidence.append({
            'session': scope.session, 'scope_start_monotonic': scope.start, 'scope_end_monotonic': scope.stop,
            'first_coordinator_monotonic': times[0], 'last_coordinator_monotonic': times[-1],
            'coordinator_count': len(rows), 'matrix_count': len(pool),
            'maximum_coordinator_gap_seconds': max(right - left for left, right in zip(times, times[1:])),
            'maximum_sampled_scope_memory_bytes': max(row['scope_memory_bytes'] for row in pool),
        })
        selected.extend(rows)
    return selected, evidence


def elapsed_seconds(value):
    seconds = 0
    for piece in value.split(':'):
        seconds = seconds * 60 + float(piece)
    return seconds


def verify_attempt(directory, phase, file, lines, scopes):
    started = int(file.stem.removeprefix(phase + '-')) / 1_000_000_000
    elapsed = [line.rsplit(': ', 1)[1] for line in lines if 'Elapsed (wall clock) time' in line]
    if len(elapsed) != 1:
        raise ValueError(f'no complete measured phase interval: {directory / file.name}')
    duration = elapsed_seconds(elapsed[0])
    if not math.isfinite(duration) or duration < 0:
        raise ValueError(f'invalid phase elapsed time: {directory / file.name}')
    ended = started + duration
    if not any(scope.unix_start <= started <= ended <= scope.unix_stop for scope in scopes):
        raise ValueError(f'measured phase did not run inside one evidenced scope: {directory / file.name}')


class InterruptedResources:
    """Validate the pinned interrupted Step 10 runtime and genuine observations."""

    RUNTIME = '20533e261485622436e5675bbb45250475d28eb6'
    SOURCE = {
        'cpst/06_tsg_module/06.4_round4/execution/matrix_measured.py':
            'f9e2c5179795236e1d87ba218ea3dd5383c5bba17d714863fb44660f5ae3170a',
        'cpst/10_cpst_dose/execution/dose_matrix.py':
            '88dbf66b2741478944940b429f3dd43d0ee52de7fc3ca8bff22943b93c7368cf',
    }
    LAUNCHER = '56d2bfc4812fed384074e4168f88cdc92cedcb6bd8dbb2a048d017d4f2210325'

    def __init__(self, bundle, pin):
        from hashlib import sha256
        from subprocess import check_output

        repo = Path(__file__).resolve().parents[3]
        root = bundle / 'matrix-audit/interruption'
        record = json.loads(pinned_file(root, 'receipt.json', pin).read_text())
        if (record['schema'] != SCHEMA or record['runtime'] != self.RUNTIME
                or record['preflight_sha256'] != digest(bundle / 'matrix/preflight.json')):
            raise ValueError('interruption receipt does not bind the pinned matrix runtime')
        # The immutable original worker reads both exact cgroup limits before its
        # ResourceMonitor starts. The original scope has no kernel peak receipt.
        for relative, expected in self.SOURCE.items():
            original = check_output(['git', '-C', repo, 'show', f'{self.RUNTIME}:{relative}'])
            if sha256(original).hexdigest() != expected:
                raise ValueError('interrupted runtime cgroup guard source differs')
        pinned_file(root, 'launcher.sh', self.LAUNCHER)
        self.scopes = scope_windows(root, record['journals'])
        settings = dict(line.split('=', 1) for line in pinned_file(
            root, 'resumed-scope-settings.txt', record['resumed_scope_settings_sha256']
        ).read_text().splitlines())
        if (Path(settings['ControlGroup']).name != record['journals'][-1]['unit']
                or (int(settings['MemoryMax']), int(settings['MemorySwapMax'])) != (12 * 1024**3, 8 * 1024**3)
                or abs(int(settings['ActiveEnterTimestampMonotonic']) / 1_000_000 - self.scopes[-1].start) > 0.001):
            raise ValueError('resumed scope cgroup settings disagree with the systemd journal')
        self.supervisor = []
        if len(record['coordinator']) != 2:
            raise ValueError('interruption proof requires two raw coordinator logs')
        for index, item in enumerate(record['coordinator']):
            rows = json_lines(pinned_file(root, item['path'], item['sha256']))
            scope = self.scopes[index]
            if not rows or not scope.start - 20 <= rows[0]['monotonic_seconds'] <= scope.start:
                raise ValueError('coordinator log does not start before its evidenced scope')
            self.supervisor.extend(rows)
        if any(right['monotonic_seconds'] <= left['monotonic_seconds']
               for left, right in zip(self.supervisor, self.supervisor[1:])):
            raise ValueError('raw coordinator logs overlap or contain nonmonotonic observations')
        self.pin = pin

    def coverage(self, matrix, final_measurement):
        selected, windows = active_samples(self.supervisor, matrix, self.scopes, final_measurement)
        limit = 12 * 1024**3
        if any(not 0 < row['scope_memory_bytes'] <= limit for row in matrix):
            raise ValueError('sampled scope memory exceeds the pinned cgroup envelope')
        return selected, {
            'interruption_receipt_sha256': self.pin,
            'samples': len(selected), 'sample_interval_seconds': 10,
            'maximum_gap_seconds': max(window['maximum_coordinator_gap_seconds'] for window in windows),
            'gaps_over_12_seconds': 0,
            'idle_window_monotonic': [self.scopes[0].stop, self.scopes[1].start],
            'active_windows': windows,
            'kernel_peak_basis': 'completed-scope-only',
            'interrupted_scope_cap_basis': 'pinned launcher and runtime cgroup guard, verified before observations',
        }

    def check_phase(self, directory, phase, file, lines):
        verify_attempt(directory, phase, file, lines, self.scopes)
