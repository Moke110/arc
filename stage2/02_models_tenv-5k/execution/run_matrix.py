"""Complete registered matrix; pure planning, strict resume and no-training acceptance."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

STEP=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(STEP/'output/api'))
METHODS={'prnet':'PRnet','transigen':'TranSiGen','mlp':'MLP_Morgan'}


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def plan():
    return [dict(adapter=a,scenario=s,fold=f,seed=t) for s in ('cold-drug','warm-start')
            for a in METHODS for f in range(5) for t in (17,42,2026)]


def require_complete(*,plan,found):
    expected={json.dumps(r,sort_keys=True) for r in plan}
    actual=[json.dumps(r,sort_keys=True) for r in found]
    if len(actual)!=len(expected) or set(actual)!=expected:raise ValueError('incomplete/duplicate matrix')


def run_path(root,run):
    return Path(root)/'runs'/run['scenario']/run['adapter']/str(run['fold'])/str(run['seed'])


def verify_manifest(root,input_digest):
    root=Path(root);manifest=json.loads((root/'manifest.json').read_text())
    if manifest['input_manifest_sha256']!=input_digest:raise ValueError('input identity changed')
    if not manifest['outputs']:raise ValueError('empty manifest')
    for name,expected in manifest['outputs'].items():
        path=root/name
        if Path(name).is_absolute() or '..' in Path(name).parts or not path.is_file() or digest(path)!=expected:
            raise ValueError('artifact digest mismatch: '+name)
    return manifest


def verify_run(root,run,input_digest,device=None):
    root=Path(root);terminal=json.loads((root/'terminal.json').read_text())
    if terminal['status']!='complete' or terminal.get('smoke') is not False or terminal.get('replay')!='exact native sampled prediction':
        raise ValueError('not a complete formal replayed run')
    config=json.loads((root/'config.json').read_text())
    for key,value in run.items():
        if config[key]!=value:raise ValueError('run identity changed')
    if config['smoke_rows'] is not None or config['max_epochs']!=1000 or config['patience']!=50 or not config['strict_lower']:
        raise ValueError('unregistered stopping/configuration')
    if device is not None and config['device']!=device:raise ValueError('device changed')
    required={'checkpoint.pt','prediction.npy','identity.csv','per_well.csv','per_drug.csv','learning_curve.csv','config.json','environment.json','terminal.json'}
    manifest=verify_manifest(root,input_digest)
    if not required<=set(manifest['outputs']):raise ValueError('missing learner artifacts')
    return config


def frozen_plan(args):
    code={str(p.relative_to(STEP)):digest(p) for directory in (STEP/'execution',STEP/'output/api')
          for p in sorted(directory.rglob('*.py')) if 'tests' not in p.parts and '__pycache__' not in p.parts}
    return dict(schema='tenv5k-matrix-v1',runs=plan(),device=args.device,input_manifest_sha256=digest(args.inputs/'manifest.json'),
        means_manifest_sha256=digest(args.means/'manifest.json'),code_sha256=code,
        sources={a:str(p.resolve()) for a,p in source_paths(args).items()},max_epochs=1000,patience=50)


def source_paths(args):
    return dict(prnet=args.prnet_source,transigen=args.transigen_source,mlp=STEP/'output/api/tenv5k')


def execute(args):
    args.candidate.mkdir(parents=True,exist_ok=True)
    frozen=frozen_plan(args);path=args.candidate/'plan.json'
    if path.exists():
        if not args.resume or json.loads(path.read_text())!=frozen:raise ValueError('resume plan/config/code identity changed')
    else:path.write_text(json.dumps(frozen,indent=2)+'\n')
    completed=0
    for run in plan():
        target=run_path(args.candidate,run)
        if target.exists():
            if not args.resume:raise ValueError('candidate run already exists')
            verify_run(target,run,frozen['input_manifest_sha256'],args.device)
            continue
        target.parent.mkdir(parents=True,exist_ok=True)
        timing=args.candidate/'resources'/('/'.join(map(str,run.values())).replace('/','-')+'.time')
        timing.parent.mkdir(exist_ok=True)
        command=[str(STEP/'execution/capped.sh'),'/usr/bin/time','-v','-o',str(timing),args.python,str(STEP/'execution'/ (run['adapter']+'.py')),
                 '--inputs',str(args.inputs),'--candidate',str(target),'--source',str(source_paths(args)[run['adapter']]),
                 '--scenario',run['scenario'],'--fold',str(run['fold']),'--seed',str(run['seed']),'--device',args.device]
        print(json.dumps(dict(start=run,command=command)),flush=True)
        with timing.with_suffix('.stdout').open('w') as stdout, timing.with_suffix('.stderr').open('w') as stderr:
            subprocess.run(command,check=True,stdout=stdout,stderr=stderr)
        verify_run(target,run,frozen['input_manifest_sha256'],args.device)
        print(json.dumps(dict(complete=run)),flush=True)
        completed+=1
        if args.limit_runs and completed>=args.limit_runs:return
    accept(args)


def evaluate(args,destination):
    """OOF first, then parent vectors, then seed metrics, then drug pairing."""
    import numpy as np
    import pandas as pd
    from tenv5k.inputs import Inputs
    from tenv5k.protocol import score_prediction
    from tenv5k.comparison import assemble_oof,seed_average,paired_comparisons,select_winner,METRICS
    data=Inputs(args.inputs);destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    input_digest=digest(args.inputs/'manifest.json')
    verify_manifest(args.means,input_digest)
    statuses=pd.read_csv(args.means/'status.csv')
    expected={(s,f,m):('not-applicable' if s=='cold-drug' and m=='Mean_drug' else 'complete') for s in ('cold-drug','warm-start') for f in range(5) for m in ('Mean','Mean_cell','Mean_drug')}
    actual={(r.scenario,int(r.fold),r.method):r.status for r in statuses.itertuples()}
    if len(statuses)!=30 or actual!=expected:raise ValueError('invalid deterministic baseline matrix')
    parent_rows=[];well_rows=[];curves=[];found=[]
    for scenario in ('cold-drug','warm-start'):
        for adapter,method in list(METHODS.items())+[(m,m) for m in ('Mean','Mean_cell','Mean_drug')]:
            if scenario=='cold-drug' and method=='Mean_drug':continue
            for seed in (17,42,2026) if adapter in METHODS else (-1,):
                folds=[]
                for fold in range(5):
                    run=dict(adapter=adapter,scenario=scenario,fold=fold,seed=seed)
                    root=run_path(args.candidate,run) if adapter in METHODS else args.means/scenario/str(fold)/method
                    if adapter in METHODS:
                        config=verify_run(root,run,input_digest,args.device);found.append(run)
                        curve=pd.read_csv(root/'learning_curve.csv')
                        if not len(curve) or len(curve)>1000 or not np.isfinite(curve[['train_loss','validation_mse','best_mse','lr']].to_numpy()).all():raise ValueError('invalid learning history')
                        terminal=json.loads((root/'terminal.json').read_text())
                        if terminal['best_epoch']!=int(curve.loc[curve.validation_mse.idxmin(),'epoch']):raise ValueError('best checkpoint history mismatch')
                        for key,value in run.items():curve[key]=value
                        curve['method']=method;curve['best_epoch']=terminal['best_epoch'];curves.append(curve)
                    identity=pd.read_csv(root/'identity.csv').profile_id.astype(str).tolist()
                    rows=data.rows(scenario,fold,'test')
                    if identity!=list(data.obs.index[rows]):raise ValueError('test order changed')
                    prediction=np.load(root/'prediction.npy',mmap_mode='r')
                    if prediction.shape!=(len(rows),5000):raise ValueError('prediction shape mismatch')
                    folds.append((identity,prediction))
                prediction=assemble_oof(data.obs.index,folds)
                wells,parents=score_prediction(data.obs,data.expression[:len(data.obs)],prediction,data.control[:len(data.obs)],prediction_ids=data.obs.index)
                for table in (wells,parents):table['scenario']=scenario;table['method']=method;table['seed']=seed
                parent_rows.append(parents);well_rows.append(wells)
    require_complete(plan=plan(),found=found)
    raw=pd.concat(parent_rows,ignore_index=True);parents=seed_average(raw)
    pairs=paired_comparisons(parents,seed=data.config['statistical_seed'])
    summaries=[]
    for (scenario,method),rows in parents.groupby(['scenario','method']):
        for metric in METRICS:
            valid=rows[metric].notna()
            summaries.append(dict(scenario=scenario,method=method,metric=metric,value=float(rows.loc[valid,metric].mean()),
                n_parents=int(valid.sum()),n_total=len(rows),n_undefined=int((~valid).sum()),zero=0.))
    raw.to_csv(destination/'per_parent_seed.csv',index=False);parents.to_csv(destination/'per_parent.csv',index=False)
    pd.concat(well_rows,ignore_index=True).to_csv(destination/'per_well.csv',index=False)
    pd.DataFrame(summaries).to_csv(destination/'metrics_values.csv',index=False)
    pairs.to_csv(destination/'paired_differences_values.csv',index=False)
    pd.concat(curves,ignore_index=True).to_csv(destination/'learning_curve_values.csv',index=False)
    parents[parents.constant_truth].to_csv(destination/'constant_truth.csv',index=False)
    (destination/'selection.json').write_text(json.dumps(select_winner(pairs),indent=2)+'\n')
    statuses.to_csv(destination/'baseline_status.csv',index=False)
    return dict(n_runs=len(found),n_wells=len(data.obs),n_parents=len(data.parents),selection=select_winner(pairs))


def accept(args):
    frozen=json.loads((args.candidate/'plan.json').read_text())
    if frozen!=frozen_plan(args):raise ValueError('acceptance plan identity changed')
    report=args.candidate/'comparison'
    result=evaluate(args,report)
    rendered=args.candidate/'figures';rendered.mkdir(exist_ok=True)
    subprocess.run([args.python,str(STEP/'figure/script/render_figures.py'),'--figure-data',str(report),'--out-dir',str(rendered)],check=True)
    second=args.candidate/'comparison-replay';evaluate(args,second)
    if {p.name:digest(p) for p in report.iterdir()}!={p.name:digest(p) for p in second.iterdir()}:raise ValueError('metric replay differs')
    redraw=args.candidate/'figures-replay'
    subprocess.run([args.python,str(STEP/'figure/script/render_figures.py'),'--figure-data',str(second),'--out-dir',str(redraw)],check=True)
    if {p.name:digest(p) for p in rendered.iterdir()}!={p.name:digest(p) for p in redraw.iterdir()}:raise ValueError('CSV figure replay differs')
    result.update(schema='tenv5k-matrix-acceptance-v1',input_manifest_sha256=digest(args.inputs/'manifest.json'),
                  plan_sha256=digest(args.candidate/'plan.json'),replay='exact CSV metrics and PNG',runs=90,
                  files={str(p.relative_to(args.candidate)):digest(p) for directory in (report,rendered) for p in sorted(directory.iterdir())})
    (args.candidate/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    if args.publish:
        publish(args)
    print(json.dumps(result),flush=True)


def publish(args):
    """Only a fully accepted candidate can replace stable products."""
    if not (args.candidate/'acceptance.json').is_file():raise ValueError('candidate not accepted')
    # Preserve public API ownership; publish independent stable payload directories.
    for source,target in [(args.candidate/'runs',STEP/'output/runs'),(args.candidate/'comparison',STEP/'output/comparison')]:
        if target.exists():
            old={str(p.relative_to(target)):digest(p) for p in target.rglob('*') if p.is_file()}
            new={str(p.relative_to(source)):digest(p) for p in source.rglob('*') if p.is_file()}
            if old!=new:raise ValueError('refuse to overwrite different accepted products')
        else:shutil.copytree(source,target)
    (STEP/'figure/data').mkdir(parents=True,exist_ok=True)
    for name in ('metrics_values.csv','paired_differences_values.csv','learning_curve_values.csv'):
        shutil.copy2(args.candidate/'comparison'/name,STEP/'figure/data'/name)
    for p in (args.candidate/'figures').iterdir():shutil.copy2(p,STEP/'figure'/p.name)
    (STEP/'log').mkdir(exist_ok=True)
    shutil.copy2(args.candidate/'acceptance.json',STEP/'log/231-acceptance.json')
    shutil.copy2(args.candidate/'plan.json',STEP/'log/231-plan.json')
    shutil.copytree(args.candidate/'resources',STEP/'log/resources',dirs_exist_ok=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan-only',action='store_true');mode.add_argument('--all',action='store_true');mode.add_argument('--accept-only',action='store_true')
    parser.add_argument('--inputs',type=Path);parser.add_argument('--means',type=Path)
    parser.add_argument('--candidate',type=Path,default=STEP/'execution/checkpoints/matrix-231')
    parser.add_argument('--prnet-source',type=Path);parser.add_argument('--transigen-source',type=Path)
    parser.add_argument('--device',default='cuda');parser.add_argument('--python',default=sys.executable)
    parser.add_argument('--resume',action='store_true');parser.add_argument('--publish',action='store_true')
    parser.add_argument('--limit-runs',type=int,help='profile this many complete formal identities, then resume unchanged')
    args=parser.parse_args()
    if args.plan_only:print(json.dumps(plan(),indent=2));return
    if any(p is None for p in (args.inputs,args.means,args.prnet_source,args.transigen_source)):parser.error('explicit inputs, means and both source paths required')
    if args.all:execute(args)
    else:accept(args)

if __name__=='__main__':main()
