"""Registered147-fit matrix with validation-first choice and strict exact replay."""
import argparse
import json
import os
from pathlib import Path
import secrets
import shutil
import subprocess
import sys
import time
import numpy as np
import pandas as pd
STEP=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(STEP/'output/api'))
from tsg_dose.protocol import METHODS,SEEDS,parent_metrics,score_prediction
from tsg_dose.comparison import choose_method,seed_average,paired_comparisons,assemble_oof
from tsg_dose.inputs import sha256,Inputs
SCENARIOS=('cold-drug','held-dose','final')
SCHEMA='tsg-dose-matrix-v1'


def read(path):return json.loads(Path(path).read_text())


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')


def plan(scenario=None):
    if scenario is not None and scenario not in SCENARIOS:raise ValueError('unknown matrix scenario')
    return [dict(scenario=s,fold=f,method=m,seed=t) for s in SCENARIOS if scenario is None or s==scenario for f in range(1 if s=='final' else 3) for m in METHODS for t in SEEDS]


def require_complete(*,plan,found):
    expected={json.dumps(r,sort_keys=True) for r in plan};actual=[json.dumps(r,sort_keys=True) for r in found]
    if len(actual)!=len(expected) or set(actual)!=expected:raise ValueError('incomplete/duplicate matrix')


def run_path(root,run):return Path(root)/'runs'/run['scenario']/run['method']/str(run['fold'])/str(run['seed'])


def registration_path(root,run):return Path(root)/'registrations'/run['scenario']/run['method']/str(run['fold'])/(str(run['seed'])+'.json')


def selection_path(root,scenario,fold):return Path(root)/'selections'/scenario/str(fold)/'receipt.json'


def validation_table(root,scenario,fold,input_digest):
    rows=[];origins={}
    for run in plan(scenario):
        if run['fold']!=fold:continue
        path=run_path(root,run);frozen=read(path/'validation_manifest.json')
        if frozen.get('schema')!='tsg-dose-validation-v1' or frozen.get('input_manifest_sha256')!=input_digest:raise ValueError('validation origin input/schema mismatch')
        for name,digest in frozen['outputs'].items():
            if Path(name).name!=name or sha256(path/name)!=digest:raise ValueError('validation artifact changed')
        config=read(path/'config.json');terminal=read(path/'terminal.json')
        if any(config[k]!=v for k,v in run.items()) or config.get('smoke_rows') is not None:raise ValueError('validation run identity mismatch')
        if terminal['status'] not in ['validation-complete','test-awaiting-replay','complete']:raise ValueError('validation not replayed')
        if frozen.get('parameter_count')!=terminal.get('parameter_count') or frozen.get('parameter_count') is None:raise ValueError('validation parameter count changed')
        scores=read(path/'validation_metrics.json')
        rows.append(dict(run,parameter_count=frozen['parameter_count'],**scores))
        origins[str(path.resolve())]=sha256(path/'validation_manifest.json')
    return pd.DataFrame(rows),origins


def lock_selection(root,scenario,fold,input_digest):
    """Seal21 replayed validation scores before the first outer-label access."""
    rows,origins=validation_table(root,scenario,fold,input_digest)
    method='no-dose' if scenario=='held-dose' else choose_method(rows)
    if scenario=='held-dose':choose_method(rows)  # validate complete seed/parameter table, without choosing an auxiliary winner
    receipt=dict(schema='tsg-dose-selection-v1',scenario=scenario,fold=fold,input_manifest_sha256=input_digest,method=method,validation_runs=origins,rule=('fixed comparison of all methods; no auxiliary method selection' if scenario=='held-dose' else 'mean three seed parent-equal validation Delta PCC; exact ties fewer parameters then frozen method order'))
    path=selection_path(root,scenario,fold)
    if path.exists():
        if read(path)!=receipt:raise ValueError('locked selection changed')
    else:
        write(path,receipt);rows.to_csv(path.parent/'validation_table.csv',index=False)
    return path


def frozen_plan(args):
    """Freeze all scientific/API/entrypoint identities once, before any fit."""
    from tsg_dose.learner import runtime_identity,environment
    from tsg_dose import model
    data=Inputs(args.inputs)
    sources=[STEP/'execution'/name for name in ['run_matrix.py','transigen.py','prepare.py','capped.sh']]+list((STEP/'output/api/tsg_dose').glob('*.py'))
    code={str(p.relative_to(STEP)):sha256(p) for p in sorted(sources)}
    transforms={}
    for scenario in SCENARIOS:
        for fold in range(1 if scenario=='final' else 3):
            doses=data.doses(data.rows(scenario,fold,'train'))
            for method in METHODS:transforms[f'{scenario}/{fold}/{method}']=model.configuration(method,model.fit_transform(method,doses))
    return dict(schema=SCHEMA,runs=plan(),device=args.device,input_manifest_sha256=sha256(args.inputs/'manifest.json'),runtime_identity=runtime_identity(),execution_code_sha256=code,source=str(args.source.resolve()),native_source_sha256=model.SOURCE_DIGEST,environment=environment(args.device),configurations=transforms,max_epochs=1000,patience=50,validation_seed=22617)


def execution_registration(run,frozen,nonce):
    return dict(run=run,device=frozen['device'],attempt_nonce=nonce,plan_sha256=__import__('hashlib').sha256(json.dumps(frozen,sort_keys=True,separators=(',',':')).encode()).hexdigest(),input_manifest_sha256=frozen['input_manifest_sha256'],origin='formal fresh training')


def verify_run(root,run,frozen,registration,*,require_test=False):
    """Verify the training origin, full digest coverage and sampled replay proof."""
    from tsg_dose.learner import verify_candidate
    from tenv5k.attempt_identity import attempt_binding
    path=run_path(root,run)
    if registration!=execution_registration(run,frozen,registration.get('attempt_nonce')):raise ValueError('execution registration changed')
    bound=attempt_binding(registration);manifest=verify_candidate(path)
    if manifest['input_manifest_sha256']!=frozen['input_manifest_sha256'] or manifest.get('execution_binding')!=bound:raise ValueError('input/attempt manifest mismatch')
    config=read(path/'config.json');terminal=read(path/'terminal.json')
    if read(path/'registration.json')!=registration or config.get('execution_binding')!=bound:raise ValueError('training attempt mismatch')
    for key,value in run.items():
        if config.get(key)!=value:raise ValueError('run identity changed')
    if config.get('runtime_identity')!=frozen['runtime_identity'] or config.get('device')!=frozen['device'] or config.get('native_source_sha256')!=frozen['native_source_sha256']:raise ValueError('runtime/source identity changed')
    expected=frozen['configurations'][f"{run['scenario']}/{run['fold']}/{run['method']}"]
    if any(config.get(k)!=v for k,v in expected.items()):raise ValueError('scientific configuration changed')
    if config.get('smoke_rows') is not None or config.get('max_epochs')!=1000 or config.get('patience')!=50 or config.get('validation_seed')!=22617:raise ValueError('unregistered fitting protocol')
    if terminal.get('smoke') is not False or terminal.get('weights_updated') is not True or terminal.get('replay')!='exact fresh-process native sampled prediction':raise ValueError('not a formal replayed fit')
    required=['validation']+(['test'] if require_test else [])
    if terminal['status'] not in (['complete'] if require_test else ['validation-complete','complete']):raise ValueError('incomplete replayed phase')
    for phase in required:
        proof=terminal.get(phase+'_replay',{})
        if not proof.get('process_token') or proof.get('pid')==terminal.get('training_pid') or proof.get('process_token')==terminal.get('training_process_token') or proof.get('attempt_nonce')!=bound['attempt_nonce']:raise ValueError('missing fresh-process replay proof')
        prediction=np.load(path/(phase+'_prediction.npy'),mmap_mode='r')
        identities=pd.read_csv(path/(phase+'_identity.csv'))
        if prediction.dtype!=np.float32 or len(prediction)!=len(identities) or identities.profile_id.duplicated().any() or not np.isfinite(prediction).all():raise ValueError('invalid serialized prediction identity/precision')
    if read(path/'environment.json')!=frozen['environment']:raise ValueError('execution environment changed')
    import torch
    checkpoint=torch.load(path/'checkpoint.pt',map_location='cpu',weights_only=False)
    if checkpoint['config']!=config or checkpoint['input_manifest_sha256']!=frozen['input_manifest_sha256']:raise ValueError('checkpoint origin mismatch')
    if any(v.is_floating_point() and v.dtype!=torch.float32 for v in checkpoint['state_dict'].values()):raise ValueError('checkpoint precision changed')
    return terminal


def event(root,kind,**data):
    path=Path(root)/'test-access-sequence.jsonl';path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a') as stream:stream.write(json.dumps(dict(event=kind,time_ns=time.time_ns(),**data),sort_keys=True)+'\n')


def launch(args,run,phase):
    path=run_path(args.candidate,run);registration=registration_path(args.candidate,run)
    key=f"{run['scenario']}-{run['fold']}-{run['method']}-{run['seed']}-{phase}"
    timing=args.candidate/'resources'/key;timing.parent.mkdir(parents=True,exist_ok=True)
    command=[str(STEP/'execution/capped.sh'),'/usr/bin/time','-v','-o',str(timing.with_suffix('.time')),args.python,str(STEP/'execution/transigen.py'),'--inputs',str(args.inputs),'--candidate',str(path),'--source',str(args.source)]
    if phase=='fit':command += ['--method',run['method'],'--scenario',run['scenario'],'--fold',str(run['fold']),'--seed',str(run['seed']),'--device',args.device,'--registration',str(registration)]
    elif phase=='replay':command += ['--replay-only']
    elif phase=='evaluate':command += ['--evaluate-only','--selection',str(selection_path(args.candidate,run['scenario'],run['fold']))]
    else:raise ValueError('unknown phase')
    print(json.dumps(dict(start=run,phase=phase)),flush=True)
    with timing.with_suffix('.stdout').open('a') as stdout,timing.with_suffix('.stderr').open('a') as stderr:subprocess.run(command,stdout=stdout,stderr=stderr,check=True)
    print(json.dumps(dict(complete=run,phase=phase)),flush=True)


def train_or_resume(args,run,frozen):
    path=run_path(args.candidate,run);registration_file=registration_path(args.candidate,run)
    if path.exists():
        if not args.resume:raise ValueError('existing origin requires --resume')
        registration=read(registration_file)
        if registration!=execution_registration(run,frozen,registration.get('attempt_nonce')):raise ValueError('resumed registration changed')
        terminal=read(path/'terminal.json') if (path/'terminal.json').exists() else {}
        if terminal.get('status')=='failed':raise ValueError('failed scientific attempt blocks release; inspect preserved evidence')
        if terminal.get('status') in ['validation-trained-awaiting-replay','test-awaiting-replay','validation-complete','complete']:
            launch(args,run,'replay');verify_run(args.candidate,run,frozen,registration,require_test=(terminal.get('status') in ['test-awaiting-replay','complete']))
            return
        # A killed process without a terminal outcome is an interrupted attempt,
        # retained verbatim before rerunning exactly the registered scientific slot.
        archive=args.candidate/'attempts'/run['scenario']/run['method']/str(run['fold'])/str(run['seed'])/registration['attempt_nonce']
        archive.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(path),str(archive))
        shutil.copy2(registration_file,archive/'driver-registration.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    registration=execution_registration(run,frozen,secrets.token_hex(32));write(registration_file,registration)
    launch(args,run,'fit');verify_run(args.candidate,run,frozen,registration)


def execute(args):
    args.candidate.mkdir(parents=True,exist_ok=True)
    frozen=frozen_plan(args);plan_path=args.candidate/'plan.json'
    if plan_path.exists():
        if not args.resume or read(plan_path)!=frozen:raise ValueError('resume input/code/configuration identity changed')
    else:write(plan_path,frozen)
    completed=0
    for scenario in ([args.scenario] if args.scenario else SCENARIOS):
        for fold in range(1 if scenario=='final' else 3):
            runs=[r for r in plan(scenario) if r['fold']==fold]
            for run in runs:
                train_or_resume(args,run,frozen);completed+=1
                if args.limit_runs and completed>=args.limit_runs:return
            receipt=lock_selection(args.candidate,scenario,fold,frozen['input_manifest_sha256'])
            event(args.candidate,'choice-locked',scenario=scenario,fold=fold,receipt_sha256=sha256(receipt))
            if scenario!='final':
                for run in runs:
                    path=run_path(args.candidate,run)
                    if read(path/'terminal.json')['status']!='complete':
                        event(args.candidate,'test-access-start',run=run,receipt_sha256=sha256(receipt));launch(args,run,'evaluate')
                    verify_run(args.candidate,run,frozen,read(registration_path(args.candidate,run)),require_test=True)
    accept(args)


def verify_test_sequence(root,scenario,fold,runs):
    receipt=selection_path(root,scenario,fold);digest=sha256(receipt)
    sequence=[json.loads(line) for line in (Path(root)/'test-access-sequence.jsonl').read_text().splitlines()]
    locked=False;seen=set()
    for entry in sequence:
        if entry['event']=='choice-locked' and entry.get('scenario')==scenario and entry.get('fold')==fold:
            if entry['receipt_sha256']!=digest:raise ValueError('choice event changed')
            locked=True
        if entry['event']=='test-access-start' and entry['run']['scenario']==scenario and entry['run']['fold']==fold:
            if not locked or entry['receipt_sha256']!=digest:raise ValueError('test access before fixed selection')
            seen.add(json.dumps(entry['run'],sort_keys=True))
    if seen!={json.dumps(r,sort_keys=True) for r in runs}:raise ValueError('missing outer test access sequence')
    for run in runs:
        stored=read(run_path(root,run)/'selection_receipt.json')
        if stored!=dict(receipt=read(receipt),sha256=digest):raise ValueError('test prediction choice binding mismatch')


def report(args,scenario,destination):
    """Complete OOF identity coverage and metric tables; one array at a time."""
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    data=Inputs(args.inputs);records=[];curves=[]
    for run in plan(scenario):
        path=run_path(args.candidate,run);terminal=read(path/'terminal.json')
        records.append(dict(run,**{k:terminal[k] for k in ['best_epoch','epochs','validation_mse','parameter_count','elapsed_seconds','peak_rss_kib']}))
        curve=pd.read_csv(path/'learning_curve.csv')
        for k,v in run.items():curve[k]=v
        curves.append(curve)
    pd.DataFrame(records).to_csv(destination/'fit_resources.csv',index=False)
    pd.concat(curves,ignore_index=True).to_csv(destination/'learning_curves.csv',index=False)
    if scenario=='final':
        receipt=read(selection_path(args.candidate,scenario,0))
        write(destination/'deployment_selection.json',dict(method=receipt['method'],seeds=list(SEEDS),role='internal-validation-selected restored checkpoints; no new independent test',checkpoints=[str(run_path(args.candidate,dict(scenario='final',fold=0,method=receipt['method'],seed=s))/'checkpoint.pt') for s in SEEDS],receipt_sha256=sha256(selection_path(args.candidate,scenario,0))))
        return
    parent_rows=[];dose_rows=[];well_rows=[];coverage=[]
    methods=list(METHODS)+(['selected'] if scenario=='cold-drug' else [])
    for method in methods:
        for seed in SEEDS:
            folds=[]
            for fold in range(3):
                actual=read(selection_path(args.candidate,scenario,fold))['method'] if method=='selected' else method
                path=run_path(args.candidate,dict(scenario=scenario,fold=fold,method=actual,seed=seed))
                identities=pd.read_csv(path/'test_identity.csv')
                expected=data.obs.index[data.rows(scenario,fold,'test')].astype(str).tolist()
                if identities.profile_id.astype(str).tolist()!=expected:raise ValueError('outer identity order differs from frozen input')
                folds.append((identities.profile_id,np.load(path/'test_prediction.npy',mmap_mode='r')))
            oof=assemble_oof(data.obs.index,folds)
            if oof.dtype!=np.float32:raise ValueError('OOF precision mismatch')
            wells,groups=score_prediction(data.obs,data.expression[:len(data.obs)],oof,data.control[:len(data.obs)],prediction_ids=data.obs.index)
            parents=parent_metrics(groups)
            for frame in [wells,groups,parents]:frame['method']=method;frame['seed']=seed
            parent_rows.append(parents);dose_rows.append(groups);well_rows.append(wells)
            coverage.append(dict(method=method,seed=seed,n_wells=len(oof),n_parents=len(parents),n_parent_dose=len(groups),folds=3))
            # Candidate reconstruction evidence. Public release retains only
            # selected final checkpoints; these predictions are not public API.
            out=destination/'oof';out.mkdir(exist_ok=True);np.save(out/f'{method}-{seed}.npy',oof)
            del oof
    parents=pd.concat(parent_rows,ignore_index=True);groups=pd.concat(dose_rows,ignore_index=True)
    wells=pd.concat(well_rows,ignore_index=True);averaged=seed_average(parents)
    parents.to_csv(destination/'per_parent_seed.csv',index=False);groups.to_csv(destination/'per_parent_dose_seed.csv',index=False)
    wells.to_csv(destination/'per_well_seed.csv',index=False);averaged.to_csv(destination/'per_parent.csv',index=False)
    pd.DataFrame(coverage).to_csv(destination/'oof_coverage.csv',index=False)
    summary=[]
    for method,frame in averaged.groupby('method',sort=False):
        summary.append(dict(method=method,**{metric:float(frame[metric].mean()) for metric in ['pcc','spearman','rmse','r2']},n_parents=len(frame),undefined_pcc=int(frame.pcc.isna().sum()),undefined_r2=int(frame.r2.isna().sum())))
    pd.DataFrame(summary).to_csv(destination/'summary.csv',index=False)
    per_dose=[]
    for (method,dose),frame in groups.groupby(['method','dose_um'],sort=False):
        row=dict(method=method,dose_um=float(dose),n_parents=int(frame.drug_parent.nunique()),n_seeds=3)
        if scenario=='held-dose':row['evaluation']={.2:'lower extrapolation',1.:'interpolation',10.:'upper extrapolation'}[float(dose)]
        for metric in ['pcc','spearman','rmse','r2']:row[metric]=float(frame.groupby('seed')[metric].mean().mean())
        row['constant_truth_count']=int(frame[frame.seed==SEEDS[0]].constant_truth.sum());row['undefined_pcc_count']=int(frame[frame.seed==SEEDS[0]].pcc.isna().sum())
        per_dose.append(row)
    pd.DataFrame(per_dose).to_csv(destination/'per_dose.csv',index=False)
    if scenario=='cold-drug':
        paired_comparisons(averaged[averaged.method.isin(['selected','no-dose'])]).to_csv(destination/'primary_paired.csv',index=False)
        paired_comparisons(averaged[averaged.method!='selected']).to_csv(destination/'descriptive_paired.csv',index=False)


def accept(args):
    frozen=read(args.candidate/'plan.json')
    if frozen!=frozen_plan(args):raise ValueError('registered runtime/input/configuration changed')
    accepted=[]
    scenarios=[args.scenario] if args.scenario else SCENARIOS
    for scenario in scenarios:
        for fold in range(1 if scenario=='final' else 3):
            runs=[r for r in plan(scenario) if r['fold']==fold]
            lock_selection(args.candidate,scenario,fold,frozen['input_manifest_sha256'])
            if scenario!='final':verify_test_sequence(args.candidate,scenario,fold,runs)
            for run in runs:
                verify_run(args.candidate,run,frozen,read(registration_path(args.candidate,run)),require_test=scenario!='final');accepted.append(run)
        report(args,scenario,args.candidate/'reports'/scenario)
    require_complete(plan=[r for r in plan() if r['scenario'] in scenarios],found=accepted)
    result=dict(schema='tsg-dose-matrix-acceptance-v1',scenarios=list(scenarios),n_formal_fits=len(accepted),runs=accepted,plan_sha256=sha256(args.candidate/'plan.json'),status='candidate complete; public release requires final acceptance',bootstrap='10000 parent-paired draws seed2262026, choices fixed; no method reselection in bootstrap')
    write(args.candidate/('acceptance-'+('-'.join(scenarios))+'.json'),result)
    print(json.dumps(result),flush=True);return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['plan','run','accept'])
    parser.add_argument('--inputs',type=Path);parser.add_argument('--candidate',type=Path,default=STEP/'execution/checkpoints/matrix')
    parser.add_argument('--source',type=Path,default=STEP.parents[1]/'pert_model/TranSiGen/source')
    parser.add_argument('--device',default='cuda');parser.add_argument('--python',default=sys.executable)
    parser.add_argument('--scenario',choices=SCENARIOS);parser.add_argument('--resume',action='store_true');parser.add_argument('--limit-runs',type=int)
    args=parser.parse_args()
    if args.action=='plan':print(json.dumps(plan(args.scenario),indent=2));return
    if args.inputs is None:parser.error('--inputs is required')
    args.inputs=args.inputs.resolve();args.candidate=args.candidate.resolve();args.source=args.source.resolve()
    # Reject code-overlapping candidates before any mkdir or stale recovery.
    if any(args.candidate==p or args.candidate in p.parents or p in args.candidate.parents for p in [STEP/'output/api',STEP/'execution/run_matrix.py',STEP/'execution/transigen.py']):raise ValueError('candidate overlaps runtime code')
    if args.action=='run':execute(args)
    else:accept(args)

if __name__=='__main__':main()
