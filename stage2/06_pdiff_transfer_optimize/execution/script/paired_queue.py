#!/usr/bin/env python3
"""Write Step 04 scheduler queues for the paired Transfer comparison."""
import argparse
from pathlib import Path

import paired_reference as reference
import paired_target as target
from artifacts import SEEDS, sha256
from learning_rate import publish
from paired_log_source import source_slot

PHASES = {
    'reference': ('verify', 'predict', 'score'),
    'log-source': ('fit',),
    **{arm: ('verify', 'fit', 'predict', 'score', 'matrix-verify') for arm in target.CHECKPOINT_NUMBERS},
}
SCRIPT = Path('execution/script')


def job_cells(arm, phase, seeds):
    if phase == 'matrix-verify':
        yield f'{arm}-matrix-verify', (arm, phase, None)
    else:
        for seed in seeds:
            yield f'{arm}-{seed}-{phase}', (arm, phase, seed)


def worker(root, runtime, arm, phase, seed):
    """Return the entry point, its sealed input arguments, and its published artifact."""
    if phase == 'matrix-verify':
        return 'paired_evidence.py', [str(root), '--arm', arm], root / arm / 'matrix.json'
    cell = root / ('log' if arm == 'log-source' else arm) / f'seed-{seed}'
    baseline = root / 'reference' / f'seed-{seed}'
    args = [phase, str(cell), '--runtime', str(runtime), '--seed', str(seed)]
    if arm == 'reference':
        artifact = {'verify': 'plan.json', 'predict': 'Transfer/prediction.json', 'score': 'scores/scores.json'}
        if phase != 'verify':
            args += ['--plan-sha256', sha256(cell / 'plan.json')]
        return 'paired_reference.py', args, cell / artifact[phase]
    if arm == 'log-source':
        args += ['--reference', str(baseline), '--plan-sha256', sha256(baseline / 'plan.json')]
        return 'paired_log_source.py', args, source_slot(runtime, seed) / 'fit/fit.json'
    args += ['--arm', arm]
    if phase == 'verify':
        args += ['--reference-sha256', sha256(baseline / 'plan.json')]
        artifact = cell / 'plan.json'
    else:
        args += ['--plan-sha256', sha256(cell / 'plan.json')]
        if phase == 'fit':
            artifact = target.target_slot(runtime, arm, seed) / 'fit/fit.json'
        else:
            artifact = cell / ('prediction/prediction.json' if phase == 'predict' else 'paired.json')
    return 'paired_target.py', args, artifact


def build_queue(root, runtime, arm, phase, seeds):
    python = runtime / 'stage2/_envs/perturbdiff/bin/python'
    jobs = []
    for job_id, cell in job_cells(arm, phase, seeds):
        script, args, _ = worker(root, runtime, *cell)
        training = phase == 'fit'
        source = arm == 'log-source'
        job = {
            'id': job_id, 'cwd': str(runtime),
            'command': [str(python), str(runtime / reference.STEP / SCRIPT / script), *args],
            'kind': 'training' if training else 'gpu' if phase == 'predict' else 'cpu',
            'memory_gib': 18 if source else 10, 'disk_gib': 5 if training else 2,
            'timeout_seconds': 15000 if training else 1800,
        }
        if training:
            job['profile'] = {'host_peak_gib': 7 if source else 5, 'gpu_peak_gib': 7}
        jobs.append(job)
    return {'record_dir': 'scheduler', 'jobs': jobs}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('arm', choices=PHASES)
    parser.add_argument('phase', choices=('verify', 'fit', 'predict', 'score', 'matrix-verify'))
    parser.add_argument('--seed', choices=SEEDS, type=int, action='append')
    parser.add_argument('--runtime', type=Path, default=reference.base.REPO)
    parser.add_argument('--queue', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.phase not in PHASES[args.arm] or (args.phase == 'matrix-verify' and args.seed):
        parser.error('phase and seed selection do not match the arm')
    seeds = args.seed or SEEDS
    if len(set(seeds)) != len(seeds):
        parser.error('duplicate seed')
    queue = build_queue(args.root.resolve(), args.runtime.resolve(), args.arm, args.phase, seeds)
    print(publish(args.queue.resolve(), queue)['path'], flush=True)


if __name__ == '__main__':
    main()
