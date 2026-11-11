"""Read a lightweight progress snapshot at fixed 30-minute wall-clock slots."""
import argparse
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import subprocess
import time


def snapshot(work, plan):
    completed, current = [], None
    for run in plan['runs']:
        fit = work / 'fits' / run['id']
        if (fit / 'run.json').exists():
            report = json.loads((fit / 'run.json').read_text())
            completed.append({'id': run['id'], 'status': report['status']})
        elif fit.exists():
            current = {'id': run['id'], 'phase': 'initializing'}
            for phase in ('tahoe', 'tenv'):
                history = fit / f'{phase}_history.json'
                if history.exists():
                    report = json.loads(history.read_text())
                    current.update(phase=phase, phase_status=report['status'],
                                   last=report['history'][-1], best_update=report['best_update'])
            if (fit / 'failure.json').exists():
                current['failure'] = json.loads((fit / 'failure.json').read_text())
    memory = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                      if line.startswith('MemAvailable:')))
    state = subprocess.run(['systemctl', '--user', 'show', 'arc-step09-round1-rerun.service',
                            '--property=ActiveState', '--property=SubState', '--property=ExecMainStatus'],
                           capture_output=True, text=True, check=True)
    return {'read_at': datetime.now(timezone(timedelta(hours=8))).isoformat(),
            'total_count': len(plan['runs']), 'complete_count': len(completed),
            'completed': completed, 'current': current,
            'available_ram_gib': memory / 1024**2,
            'service': dict(line.split('=', 1) for line in state.stdout.splitlines())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--first-at', required=True, help='timezone-aware ISO timestamp')
    args = parser.parse_args()
    due = datetime.fromisoformat(args.first_at)
    if due.tzinfo is None:
        parser.error('--first-at must include timezone')
    plan = json.loads(args.plan.read_text())
    while True:
        # Sleeping does not read training progress; the OS process keeps this
        # deadline even while the agent is waiting for a tool result.
        while due.timestamp() > time.time():
            time.sleep(max(0, min(30, due.timestamp() - time.time())))
        result = snapshot(args.work, plan)
        result['scheduled_at'] = due.isoformat()
        with (args.work / 'monitor-snapshots.jsonl').open('a') as stream:
            stream.write(json.dumps(result) + '\n')
        print(json.dumps(result), flush=True)
        if result['service']['ActiveState'] != 'active':
            break
        due += timedelta(minutes=30)


if __name__ == '__main__':
    main()
