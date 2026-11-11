"""Verify paired-run scheduler receipts independently of retrieval reports."""
import importlib.util
from pathlib import Path
import paired_reference as reference
from artifacts import SEEDS, read_json, sha256
from formal_matrix import bind
from paired_queue import PHASES, SCRIPT, job_cells, worker
STEP = reference.STEP


def relative(binding, runtime):
    return {**binding, 'path': str(Path(binding['path']).relative_to(runtime))}


def load_scheduler():
    # Step 03 also provides schedule.py. Use the Step 04 worker by file identity.
    path = reference.base.REPO / reference.TRANSFER.STEP04 / SCRIPT / 'schedule.py'
    spec = importlib.util.spec_from_file_location('step06_step04_scheduler', path)
    scheduler = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scheduler)
    return scheduler


def verify_job(job, cell, root, runtime, entry, directory, number):
    script, args, artifact = worker(root, runtime, *cell)
    python = runtime / 'stage2/_envs/perturbdiff/bin/python'
    identity = entry['identity']
    command = identity['command']
    script_suffix = (STEP / SCRIPT / script).parts
    if (len(command) < 2 or identity['cwd'] != str(runtime) or command[0] != str(python)
            or not Path(command[1]).is_absolute()
            or Path(command[1]).parts[-len(script_suffix):] != script_suffix
            or command[2:] != args):
        raise ValueError(f'{job} queued command differs from the expected cell or sealed plan')
    attempts = entry['attempts']
    if entry['status'] != 'completed' or not attempts:
        raise ValueError(f'{job} needs a completed scheduler attempt')
    attempt = attempts[-1]
    if (attempt['outcome'] != 'completed' or attempt['returncode'] != 0
            or any(not number(item['started_at'], 1) or not number(item['finished_at'], item['started_at'])
                   for item in attempts)):
        raise ValueError(f'{job} needs a successful final attempt with recorded timing')
    log = directory / f'{job}.attempt{len(attempts)}.log'
    if attempt['log'] != str(log):
        raise ValueError(f'{job} final attempt log differs from the scheduler slot')
    receipt = f'{artifact} SHA-256 {sha256(artifact)}'
    with log.open() as stream:
        if not any(line.rstrip() == receipt for line in stream):
            raise ValueError(f'{job} completed attempt did not publish the sealed artifact')
    return relative(bind(log), runtime)


def resources(root, runtime):
    scheduler = load_scheduler()
    expected = dict(job_cell for arm, phases in PHASES.items() for phase in phases
                    for job_cell in job_cells(arm, phase, SEEDS))
    files = {}
    ledger = {}
    measurements = ('memory_max_gib', 'memory_swap_max_gib', 'peak_host_rss_gib',
                    'peak_workstation_gpu_used_gib', 'peak_gpu_utilization_percent')
    for status in sorted((root / 'schedule').glob('*/scheduler/status.json')):
        jobs = read_json(status)['jobs']
        tracked = expected.keys() & jobs.keys()
        if not tracked:
            continue
        slot = status.parent.parent
        queue = slot / 'queue.json'
        record_dir, queued = scheduler.load_queue(queue)
        if record_dir != status.parent:
            raise ValueError(f'{slot.name} scheduler status is outside the queued record directory')
        jobs = scheduler.load_status(record_dir, queued)['jobs']
        files[slot.name] = [relative(bind(path), runtime) for path in
                            (queue, status, status.with_name('telemetry.jsonl'))]
        for job in sorted(tracked):
            entry = jobs[job]
            if job in ledger:
                raise ValueError(f'{job} needs one completed scheduler record')
            log = verify_job(job, expected[job], root, runtime, entry, record_dir, scheduler.number)
            ledger[job] = {'slot': slot.name, 'completed_log': log, 'attempts': [
                {'outcome': item['outcome'], 'wall_seconds': item['finished_at'] - item['started_at'],
                 **{key: item.get(key) for key in measurements}}
                for item in entry['attempts']
            ]}
    if missing := expected.keys() - ledger.keys():
        raise ValueError(f'missing scheduler evidence: {", ".join(sorted(missing))}')
    return {'jobs': ledger, 'files': files}
