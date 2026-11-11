"""Serial, capped execution of the verified ADR 0060 24-slot plan."""
import argparse
import fcntl
import json
from pathlib import Path
import subprocess
import sys

from round_contract import sha256, read_json, write_json

SCRIPT = Path(__file__).resolve().parent
MODELS = {f'XPert_{kind}{size}' for kind in ('G', 'MLP') for size in (512, 1024, 2048)}


def verify_plan(plan):
    expected = {(m, f, s) for m in MODELS for f in (0, 1) for s in (17, 42)}
    runs = plan['runs']
    if (plan['schema'] != 'step09-round1-rerun-plan-v1' or plan['fixture']
            or len(runs) != 24 or {(r['model'], r['fold'], r['seed']) for r in runs} != expected
            or len({r['id'] for r in runs}) != 24):
        raise ValueError('formal execution requires the verified complete real-data plan')
    for run in runs:
        for side in ('predictor', 'evaluator'):
            if sha256(Path(run[side]) / 'manifest.json') != run[f'{side}_sha256']:
                raise ValueError(f"{run['id']}: changed {side} manifest")


def completed(run, fit, policy):
    report = read_json(fit / 'run.json')
    if (report['status'] != 'complete' or report['model'] != run['model']
            or report['fold'] != run['fold'] or report['seed'] != run['seed']
            or report['predictor_sha256'] != run['predictor_sha256']
            or report['policy_sha256'] != sha256(policy)):
        raise ValueError(f"{run['id']}: completed fit identity differs")
    if set(p.name for p in fit.glob('*.pt')) != {'best.pt'}:
        raise ValueError('rerun must retain exactly one final checkpoint')
    if sha256(fit / 'best.pt') != report['phases']['tenv']['checkpoint_sha256']:
        raise ValueError('final checkpoint differs from selected tenv weights')
    return report


def execute(plan_path, policy, work, selected=None, smoke=False):
    plan = read_json(plan_path)
    verify_plan(plan)
    recipe = read_json(policy)
    if (recipe['kind'] == 'qualification') != smoke:
        raise ValueError('qualification and formal output must remain separate')
    work.mkdir(parents=True, exist_ok=True)
    with (work / '.queue.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for run in plan['runs']:
            if selected and run['id'] not in selected:
                continue
            fit = work / 'fits' / run['id']
            if (fit / 'run.json').exists():
                completed(run, fit, policy)
                continue
            if fit.exists():
                raise ValueError(f'incomplete fit requires explicit recovery outside queue: {fit}')
            resources = work / 'resources' / run['id']
            print(f"Starting {run['id']}", flush=True)
            args = [str(SCRIPT / 'run_capped.sh'), str(resources), sys.executable,
                    str(SCRIPT / 'xpert.py'), 'fit', '--model', run['model'],
                    '--predictor', run['predictor'], '--policy', str(policy),
                    '--fold', str(run['fold']), '--seed', str(run['seed']), '--output', str(fit)]
            result = subprocess.run(args)
            write_json(work / 'queue-status.json', {'run': run['id'], 'exit_code': result.returncode,
                       'plan_sha256': sha256(plan_path), 'smoke': smoke})
            if result.returncode:
                raise RuntimeError(f"{run['id']} failed; inspect {resources}")
            completed(run, fit, policy)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--policy', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--runs', nargs='+')
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()
    execute(args.plan, args.policy, args.work, args.runs, args.smoke)


if __name__ == '__main__':
    main()
