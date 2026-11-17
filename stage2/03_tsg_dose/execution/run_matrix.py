"""Round1: 126 validation-first cold-drug/held-dose fits; no deployment role."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

STEP=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(STEP/'output/api'))

def _frozen(name):
    spec=importlib.util.spec_from_file_location('_dose_'+name,STEP/'execution/historical'/f'{name}_v1.py')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.STEP=STEP
    return module

_historical=_frozen('run_matrix')
_active=_frozen('run_matrix')
SCENARIOS=('cold-drug','held-dose')
_active.SCENARIOS=SCENARIOS

def plan(scenario=None):
    if scenario is not None and scenario not in SCENARIOS:
        raise ValueError('round1 excludes final/deployment scenarios')
    return _active.plan(scenario)

# These scientific mechanisms retain their accepted implementation and identity.
for _name in ['require_complete','lock_selection','verify_test_sequence','write','event',
              'execution_registration','verify_run','train_or_resume','registration_path',
              'run_path','selection_path','read']:
    globals()[_name]=getattr(_active,_name)


def historical_code_mapping():
    """Resolve original execution paths to byte-identical retained source."""
    return {f'execution/{name}.py':f'execution/historical/{name}_v1.py'
            for name in ['run_matrix','transigen']}


def verify_frozen_plan(args,frozen):
    if len(frozen.get('runs',[]))==147:
        expected=_historical.frozen_plan(args)
        for original,preserved in historical_code_mapping().items():
            expected['execution_code_sha256'][original]=_active.sha256(STEP/preserved)
    else:
        expected=_active.frozen_plan(args)
    if frozen!=expected:
        raise ValueError('registered runtime/input/configuration or preserved code changed')


def accept(args):
    """Verify retained replay evidence without final fits or model publication.

    Historical scientific report tables are reused, never overwritten. A new
    round1 candidate may generate those tables after every retained origin passes.
    """
    frozen=read(args.candidate/'plan.json')
    verify_frozen_plan(args,frozen)
    scenarios=[args.scenario] if args.scenario else list(SCENARIOS)
    accepted=[]
    for scenario in scenarios:
        plan(scenario)
        for fold in range(3):
            runs=[r for r in plan(scenario) if r['fold']==fold]
            receipt=selection_path(args.candidate,scenario,fold)
            if not receipt.exists():raise ValueError('missing sealed choice receipt')
            before=_active.sha256(receipt)
            lock_selection(args.candidate,scenario,fold,frozen['input_manifest_sha256'])
            if _active.sha256(receipt)!=before:raise ValueError('sealed choice changed')
            verify_test_sequence(args.candidate,scenario,fold,runs)
            for run in runs:
                verify_run(args.candidate,run,frozen,read(registration_path(args.candidate,run)),require_test=True)
                accepted.append(run)
        if len(frozen['runs'])==126:
            _active.report(args,scenario,args.candidate/'reports'/scenario)
        elif not (args.candidate/'reports'/scenario/'summary.csv').is_file():
            raise ValueError('historical scientific report missing')
    require_complete(plan=[r for r in plan() if r['scenario'] in scenarios],found=accepted)
    result=dict(schema='tsg-dose-round1-acceptance-v1',scenarios=scenarios,
                n_formal_fits=len(accepted),runs=accepted,
                plan_sha256=_active.sha256(args.candidate/'plan.json'),
                historical_code_mapping=historical_code_mapping(),
                status='optimization evidence complete; no downstream model weights',
                replay='retained exact fresh-process validation and test evidence verified')
    write(getattr(args,'receipt',None) or args.candidate/('round1-acceptance-'+('-'.join(scenarios))+'.json'),result)
    return result

_new_frozen_plan=_active.frozen_plan
def frozen_plan(args):
    result=_new_frozen_plan(args)
    result['execution_code_sha256'].update({preserved:_active.sha256(STEP/preserved) for preserved in historical_code_mapping().values()})
    return result
_active.frozen_plan=frozen_plan
_active.accept=accept

def execute(args):
    if args.scenario is not None:plan(args.scenario)
    return _active.execute(args)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['plan','run','accept'])
    parser.add_argument('--inputs',type=Path)
    parser.add_argument('--receipt',type=Path,help='separate round1 acceptance receipt; historical records stay unchanged')
    parser.add_argument('--candidate',type=Path,default=STEP/'execution/checkpoints/round1')
    parser.add_argument('--source',type=Path,required=False)
    parser.add_argument('--device',default='cuda')
    parser.add_argument('--python',default=sys.executable)
    parser.add_argument('--scenario',choices=SCENARIOS)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--limit-runs',type=int)
    args=parser.parse_args()
    if args.action=='plan':print(json.dumps(plan(args.scenario),indent=2));return
    if args.inputs is None or args.source is None:parser.error('--inputs and --source are required')
    args.inputs=args.inputs.resolve();args.candidate=args.candidate.resolve();args.source=args.source.resolve()
    if any(args.candidate==p or args.candidate in p.parents or p in args.candidate.parents
           for p in [STEP/'output/api',STEP/'execution/run_matrix.py',STEP/'execution/transigen.py',STEP/'execution/historical']):
        raise ValueError('candidate overlaps runtime code')
    print(json.dumps(execute(args) if args.action=='run' else accept(args)))

if __name__=='__main__':main()
