"""Complete registered matrix; pure planning, strict resume and no-training acceptance."""
import argparse
import ast
import hashlib
import json
import os
import secrets
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

STEP=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(STEP/'output/api'))
sys.path.insert(0,str(STEP/'execution'))
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
    methods={}
    for adapter in METHODS:
        names=['execution/'+adapter+'.py']+['output/api/tenv5k/'+name+'.py' for name in
            ('model_'+adapter,'learner','learning','inputs','protocol','native_import','attempt_identity')]
        identity={'prnet':'prnet-native-f32-v1','transigen':'transigen-native-f32-v1','mlp':'mlp-morgan-project-v1'}[adapter]
        tree=ast.parse((STEP/'output/api/tenv5k'/('model_'+adapter+'.py')).read_text())
        for node in tree.body:
            if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='CONFIG' for t in node.targets) and isinstance(node.value,ast.Call):
                for keyword in node.value.keywords:
                    if keyword.arg=='method_identity':identity=ast.literal_eval(keyword.value)
        methods[adapter]=dict(scientific_identity=identity,code_sha256={n:code[n] for n in names})
    shared=dict(input_manifest_sha256=digest(args.inputs/'manifest.json'),
                means_manifest_sha256=digest(args.means/'manifest.json'),
                evaluation_code_sha256={n:code[n] for n in ('output/api/tenv5k/protocol.py','output/api/tenv5k/comparison.py')})
    return dict(schema='tenv5k-matrix-v2',shared=shared,methods=methods,runs=plan(),device=args.device,input_manifest_sha256=digest(args.inputs/'manifest.json'),
        means_manifest_sha256=digest(args.means/'manifest.json'),code_sha256=code,
        sources={a:str(p.resolve()) for a,p in source_paths(args).items()},max_epochs=1000,patience=50)


def import_run(args,run):
    if run not in plan() or run["adapter"]!="prnet" or run["scenario"]!="cold-drug":
        raise ValueError("only the 15 registered cold PRnet identities may be imported")
    from import_matrix import import_original
    return import_original(args,run,import_dependencies())


def import_dependencies():
    from types import SimpleNamespace
    return SimpleNamespace(plan=plan,digest=digest,verify_run=verify_run,frozen_plan=frozen_plan,
                           run_path=run_path,tree_identity=_tree_identity,step=STEP)


def source_paths(args):
    return dict(prnet=args.prnet_source,transigen=args.transigen_source,mlp=STEP/'output/api/tenv5k')


def registration_path(candidate,run):
    return Path(candidate)/"registrations"/run["scenario"]/run["adapter"]/str(run["fold"])/(str(run["seed"])+".json")


def execution_registration(run,frozen,attempt_nonce=None):
    result=dict(run=run,method=frozen["methods"][run["adapter"]],shared=frozen["shared"],
                runtime_code_sha256=frozen["code_sha256"],device=frozen["device"],sources=frozen["sources"])
    if attempt_nonce is not None:result["attempt_nonce"]=attempt_nonce
    return result


def verify_registration(registration,run,frozen):
    from tenv5k.attempt_identity import attempt_binding
    if not isinstance(registration,dict) or registration!=execution_registration(run,frozen,registration.get("attempt_nonce")):
        raise ValueError("new training lacks preregistered execution identity")
    return attempt_binding(registration)


def record_execution(root,run,frozen,registration=None):
    root=Path(root)
    if (root/'import-provenance').exists() or (root/'import-receipt.json').exists():
        raise ValueError('imported evidence cannot be relabelled as newly trained')
    config=verify_run(root,run,frozen['input_manifest_sha256'],frozen['device'])
    method=frozen['methods'][run['adapter']]
    if run['adapter']=='transigen' and config.get('method_identity','transigen-native-f32-v1')!=method['scientific_identity']:
        raise ValueError('trained numerical method identity changed')
    verify_registration(registration,run,frozen)
    from tenv5k.attempt_identity import verify_artifact_binding
    verify_artifact_binding(root,registration)
    receipt=dict(schema='tenv5k-execution-v1',origin='trained',run=run,method=method,shared=frozen['shared'],
                 registration=registration,
                 manifest_sha256=digest(root/'manifest.json'))
    (root/'execution-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


def verify_execution(root,run,frozen):
    root=Path(root)
    if (root/'import-receipt.json').exists():
        from import_matrix import verify_import
        return verify_import(root,run,frozen,import_dependencies())
    receipt=json.loads((root/'execution-receipt.json').read_text())
    registration=receipt.get('registration')
    verify_registration(registration,run,frozen)
    from tenv5k.attempt_identity import verify_artifact_binding
    verify_artifact_binding(root,registration)
    expected=dict(schema='tenv5k-execution-v1',origin='trained',run=run,method=frozen['methods'][run['adapter']],
                  shared=frozen['shared'],registration=registration,manifest_sha256=digest(root/'manifest.json'))
    if receipt!=expected:raise ValueError('trained execution identity changed')
    config=verify_run(root,run,frozen['input_manifest_sha256'],frozen['device'])
    if run['adapter']=='transigen' and config.get('method_identity','transigen-native-f32-v1')!=receipt['method']['scientific_identity']:
        raise ValueError('trained numerical method identity changed')
    return receipt


def recover_existing(args,run,input_digest,replay_only,registration=None):
    """Resume exact identities; retain interrupted bytes before a clean retry."""
    target=run_path(args.candidate,run)
    terminal_path=target/'terminal.json'
    try:terminal=json.loads(terminal_path.read_text()) if terminal_path.exists() else {}
    except json.JSONDecodeError:terminal={}
    slot=registration_path(args.candidate,run)
    if registration is None and slot.exists():registration=json.loads(slot.read_text())
    if registration is not None and terminal.get('status') in ('complete','trained-awaiting-fresh-process-replay'):
        from tenv5k.attempt_identity import verify_artifact_binding
        verify_artifact_binding(target,registration,require_manifest=terminal.get('status')=='complete')
    if terminal.get('status')=='complete':
        verify_run(target,run,input_digest,args.device)
        return 'complete'
    config_path=target/'config.json'
    try:config=json.loads(config_path.read_text()) if config_path.exists() else None
    except json.JSONDecodeError:config=None  # Interrupted write: retained verbatim below.
    if config is not None:
        expected={**run,'device':args.device,'smoke_rows':None,'max_epochs':1000,'patience':50,'strict_lower':True}
        if any(config.get(k)!=v for k,v in expected.items()):raise ValueError('cannot recover a different/smoke run identity')
    attempts=Path(args.candidate)/'attempts'/run['scenario']/run['adapter']/str(run['fold'])/str(run['seed'])
    attempts.mkdir(parents=True,exist_ok=True)
    number=1
    while (attempts/('attempt-'+str(number).zfill(4))).exists():number+=1
    attempt=attempts/('attempt-'+str(number).zfill(4))
    if terminal.get('status')=='trained-awaiting-fresh-process-replay':
        if config is None:raise ValueError('trained run lacks registered config')
        attempt.mkdir()
        if registration is not None:(attempt/'execution-registration.json').write_text(json.dumps(registration,indent=2)+'\n')
        shutil.copy2(terminal_path,attempt/'terminal.before-replay.json')
        (attempt/'recovery.json').write_text(json.dumps(dict(action='replay-only',run=run,input_manifest_sha256=input_digest),indent=2)+'\n')
        try:
            replay_only()
            verify_run(target,run,input_digest,args.device)
        except BaseException as error:
            # Preserve original evidence and allow a subsequent same-identity retry.
            if terminal_path.exists():shutil.copy2(terminal_path,attempt/'terminal.after-failed-replay.json')
            terminal_path.write_text(json.dumps(dict(status='failed',error='replay recovery: '+str(error)),indent=2)+'\n')
            raise
        return 'complete'
    # A complete-but-corrupt run never reaches this branch; partial/failed runs
    # restart under the identical frozen plan and keep their entire former tree.
    os.rename(target,attempt)
    if registration is not None:(attempt/'execution-registration.json').write_text(json.dumps(registration,indent=2)+'\n')
    (attempt/'recovery.json').write_text(json.dumps(dict(action='restart-same-identity',run=run,input_manifest_sha256=input_digest),indent=2)+'\n')
    return 'train'


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
            if (target/'import-provenance').exists() and not (target/'import-receipt.json').exists():
                raise ValueError('import receipt missing; cannot resume as newly trained')
            if (target/'import-receipt.json').exists():
                verify_execution(target,run,frozen)
                continue
            if (target/'execution-receipt.json').exists():
                verify_execution(target,run,frozen)
                continue
            registration=json.loads(registration_path(args.candidate,run).read_text())
            verify_registration(registration,run,frozen)
            replay_command=[args.python,str(STEP/'execution'/(run['adapter']+'.py')),
                            '--inputs',str(args.inputs),'--candidate',str(target),
                            '--source',str(source_paths(args)[run['adapter']]),'--replay-only',
                            '--registration',str(registration_path(args.candidate,run))]
            action=recover_existing(args,run,frozen['input_manifest_sha256'],
                                    lambda:subprocess.run([str(STEP/'execution/capped.sh'),*replay_command],check=True),registration)
            if action=='complete':
                record_execution(target,run,frozen,registration)
                continue
        target.parent.mkdir(parents=True,exist_ok=True)
        registration=execution_registration(run,frozen,secrets.token_hex(32))
        registration_file=registration_path(args.candidate,run)
        registration_file.parent.mkdir(parents=True,exist_ok=True)
        if registration_file.exists():
            history=registration_file.parent/'history'/str(run['seed'])/(digest(registration_file)+'.json')
            history.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(registration_file,history)
        registration_file.write_text(json.dumps(registration,indent=2)+'\n')
        timing=args.candidate/'resources'/('/'.join(map(str,run.values())).replace('/','-')+'.time')
        timing.parent.mkdir(exist_ok=True)
        command=[str(STEP/'execution/capped.sh'),'/usr/bin/time','-v','-o',str(timing),args.python,str(STEP/'execution'/ (run['adapter']+'.py')),
                 '--inputs',str(args.inputs),'--candidate',str(target),'--source',str(source_paths(args)[run['adapter']]),
                 '--scenario',run['scenario'],'--fold',str(run['fold']),'--seed',str(run['seed']),'--device',args.device,
                 '--registration',str(registration_file)]
        print(json.dumps(dict(start=run,command=command)),flush=True)
        with timing.with_suffix('.stdout').open('w') as stdout, timing.with_suffix('.stderr').open('w') as stderr:
            subprocess.run(command,check=True,stdout=stdout,stderr=stderr)
        verify_run(target,run,frozen['input_manifest_sha256'],args.device)
        record_execution(target,run,frozen,registration)
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
    frozen=json.loads((args.candidate/'plan.json').read_text())
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
                        verify_execution(root,run,frozen)
                        config=verify_run(root,run,input_digest,args.device);found.append(run)
                        curve=pd.read_csv(root/'learning_curve.csv')
                        if not len(curve) or len(curve)>1000 or not np.isfinite(curve[['train_loss','validation_mse','best_mse','lr']].to_numpy()).all():raise ValueError('invalid learning history')
                        terminal=json.loads((root/'terminal.json').read_text())
                        if terminal['best_epoch']!=int(curve.loc[curve.validation_mse.idxmin(),'epoch']):raise ValueError('best checkpoint history mismatch')
                        for key,value in run.items():curve[key]=value
                        curve['method_identity']=frozen['methods'][adapter]['scientific_identity']
                        curve['method']=method;curve['best_epoch']=terminal['best_epoch'];curve['status']='complete';curves.append(curve)
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
    method_labels={name:frozen['methods'][adapter]['scientific_identity'] for adapter,name in METHODS.items()}
    method_labels.update({name:'deterministic-'+name for name in ('Mean','Mean_cell','Mean_drug')})
    raw['method_identity']=raw.method.map(method_labels)
    parents['method_identity']=parents.method.map(method_labels)
    pairs=paired_comparisons(parents,seed=data.config['statistical_seed'])
    pairs['candidate_identity']=pairs.candidate.map(method_labels)
    pairs['reference_identity']=pairs.reference.map(method_labels)
    for column in ('candidate_constant_prediction','reference_constant_prediction'):
        pairs[column]=0
    pairs['constant_truth_ids']=''
    pairs['unavailable_reference']=''
    pairs['status']='complete'
    for index,pair in pairs.iterrows():
        scenario_rows=parents[parents.scenario==pair.scenario]
        pairs.at[index,'candidate_constant_prediction']=int(scenario_rows[scenario_rows.method==pair.candidate].constant_prediction_count.sum())
        pairs.at[index,'reference_constant_prediction']=int(scenario_rows[scenario_rows.method==pair.reference].constant_prediction_count.sum())
        pairs.at[index,'constant_truth_ids']=';'.join(sorted(scenario_rows[scenario_rows.constant_truth].drug_parent.unique()))
        if pair.scenario=='cold-drug':pairs.at[index,'unavailable_reference']='Mean_drug: unseen test drug'
    summaries=[]
    for (scenario,method),rows in parents.groupby(['scenario','method']):
        for metric in METRICS:
            valid=rows[metric].notna()
            summaries.append(dict(scenario=scenario,method=method,metric=metric,value=float(rows.loc[valid,metric].mean()),
                n_parents=int(valid.sum()),n_total=len(rows),n_undefined=int((~valid).sum()),zero=0.,
                applicable=True,status='complete',n_constant_truth=int(rows.constant_truth.sum()),
                n_constant_prediction=int(rows.constant_prediction_count.sum()),
                constant_truth_ids=';'.join(sorted(rows.loc[rows.constant_truth,'drug_parent'])),
                reason='truth constant: correlations/R2 undefined; error metrics retained'))
    # Inapplicability is a recorded report row, never an invisible missing method.
    cold=parents[(parents.scenario=='cold-drug')&(parents.method=='Mean')]
    for metric in METRICS:
        summaries.append(dict(scenario='cold-drug',method='Mean_drug',metric=metric,value=float('nan'),
            n_parents=0,n_total=len(cold),n_undefined=len(cold),zero=0.,applicable=False,status='not-applicable',
            n_constant_truth=int(cold.constant_truth.sum()),n_constant_prediction=0,
            constant_truth_ids=';'.join(sorted(cold.loc[cold.constant_truth,'drug_parent'])),reason='unseen test drug; no fallback'))
    raw.to_csv(destination/'per_parent_seed.csv',index=False);parents.to_csv(destination/'per_parent.csv',index=False)
    well_table=pd.concat(well_rows,ignore_index=True);well_table['method_identity']=well_table.method.map(method_labels)
    well_table.to_csv(destination/'per_well.csv',index=False)
    summary=pd.DataFrame(summaries);summary['method_identity']=summary.method.map(method_labels)
    summary.to_csv(destination/'metrics_values.csv',index=False)
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
    origins=[verify_execution(run_path(args.candidate,r),r,frozen) for r in plan()]
    result.update(execution_origins=[dict(run=r['run'],origin=r['origin']) for r in origins],
                  method_identities=frozen['methods'],shared_identity=frozen['shared'],
                  schema='tenv5k-matrix-acceptance-v2',input_manifest_sha256=digest(args.inputs/'manifest.json'),
                  plan_sha256=digest(args.candidate/'plan.json'),replay='exact CSV metrics and PNG',runs=90,
                  files={str(p.relative_to(args.candidate)):digest(p) for directory in (report,rendered) for p in sorted(directory.iterdir())})
    (args.candidate/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    if args.publish:
        publish(args)
    print(json.dumps(result),flush=True)


def _tree_identity(path):
    path=Path(path)
    if path.is_file(): return {'.':digest(path)}
    return {str(p.relative_to(path)):digest(p) for p in path.rglob('*') if p.is_file()}


def _rollback_publication(journal):
    """Restore the previous destinations, including after interrupted renames."""
    state=json.loads(journal.read_text())
    for entry in reversed(state['entries']):
        target=Path(entry['target']);backup=Path(entry['backup'])
        if backup.exists():
            if target.is_dir():shutil.rmtree(target)
            elif target.exists():target.unlink()
            target.parent.mkdir(parents=True,exist_ok=True)
            os.rename(backup,target)
        elif not entry['existed'] and not Path(entry['staged']).exists():
            if target.is_dir():shutil.rmtree(target)
            elif target.exists():target.unlink()
    journal.unlink()
    shutil.rmtree(state['transaction'],ignore_errors=True)


def publish(args):
    """Preflight all destinations; stage verified copies; rollback any install failure."""
    if not (args.candidate/'acceptance.json').is_file():raise ValueError('candidate not accepted')
    acceptance=json.loads((args.candidate/'acceptance.json').read_text())
    if acceptance.get('schema')=='tenv5k-matrix-acceptance-v2':
        frozen=json.loads((args.candidate/'plan.json').read_text())
        if acceptance['plan_sha256']!=digest(args.candidate/'plan.json'):raise ValueError('accepted plan changed')
        for run in plan():verify_execution(run_path(args.candidate,run),run,frozen)
    journal=STEP.parent/'.tenv5k-publication-journal.json'
    if journal.exists():_rollback_publication(journal)
    destinations=[(args.candidate/'runs',STEP/'output/runs',True),
                  (args.candidate/'comparison',STEP/'output/comparison',True)]
    destinations += [(args.candidate/'comparison'/name,STEP/'figure/data'/name,False)
                     for name in ('metrics_values.csv','paired_differences_values.csv','learning_curve_values.csv')]
    destinations += [(p,STEP/'figure'/p.name,False) for p in sorted((args.candidate/'figures').iterdir())]
    destinations += [(args.candidate/'acceptance.json',STEP/'log/231-acceptance.json',False),
                     (args.candidate/'plan.json',STEP/'log/231-plan.json',False),
                     (args.candidate/'resources',STEP/'log/resources',False)]
    # All refusal paths precede even the first stable-path mutation.
    for source,target,immutable in destinations:
        if not source.exists():raise ValueError('missing publication source: '+str(source))
        if target.exists() and source.is_dir()!=target.is_dir():raise ValueError('publication destination type changed: '+str(target))
        if immutable and target.exists() and _tree_identity(target)!=_tree_identity(source):
            raise ValueError('refuse to overwrite different accepted products: '+str(target))
    transaction=Path(tempfile.mkdtemp(prefix='.tenv5k-publish-',dir=STEP.parent))
    entries=[]
    try:
        for index,(source,target,_) in enumerate(destinations):
            staged=transaction/('staged-'+str(index));backup=transaction/('backup-'+str(index))
            if source.is_dir():shutil.copytree(source,staged)
            else:shutil.copy2(source,staged)
            if _tree_identity(staged)!=_tree_identity(source):raise ValueError('staged publication digest differs')
            entries.append(dict(target=str(target),staged=str(staged),backup=str(backup),existed=target.exists()))
        # Journal survives a terminated installer; the next attempt restores first.
        pending=transaction/'journal.pending'
        with pending.open('w') as stream:
            stream.write(json.dumps(dict(transaction=str(transaction),entries=entries),indent=2)+'\n')
            stream.flush();os.fsync(stream.fileno())
        os.rename(pending,journal)
        for entry in entries:
            target=Path(entry['target']);target.parent.mkdir(parents=True,exist_ok=True)
            if target.exists():os.replace(target,entry['backup'])
            os.replace(entry['staged'],target)
        journal.unlink()
    except BaseException:
        if journal.exists():
            try:
                _rollback_publication(journal)
            except BaseException as recovery_error:
                raise RuntimeError('Publication rollback failed; retained journal '+str(journal)+
                                   ' and transaction '+str(transaction)+' for the next recovery attempt') from recovery_error
        raise
    finally:
        # A remaining journal may reference the only accepted backups. Never
        # clean its transaction until installation or full rollback succeeds.
        if not journal.exists():shutil.rmtree(transaction,ignore_errors=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--import-cold-prnet',action='store_true');mode.add_argument('--plan-only',action='store_true');mode.add_argument('--all',action='store_true');mode.add_argument('--accept-only',action='store_true')
    parser.add_argument('--inputs',type=Path);parser.add_argument('--means',type=Path)
    parser.add_argument('--candidate',type=Path,default=STEP/'execution/checkpoints/matrix-231')
    parser.add_argument('--prnet-source',type=Path);parser.add_argument('--transigen-source',type=Path)
    parser.add_argument('--import-original',type=Path);parser.add_argument('--import-runtime',type=Path)
    parser.add_argument('--compatibility-audit',type=Path)
    parser.add_argument('--device',default='cuda');parser.add_argument('--python',default=sys.executable)
    parser.add_argument('--resume',action='store_true');parser.add_argument('--publish',action='store_true')
    parser.add_argument('--limit-runs',type=int,help='profile this many complete formal identities, then resume unchanged')
    args=parser.parse_args()
    if args.plan_only:print(json.dumps(plan(),indent=2));return
    if any(p is None for p in (args.inputs,args.means,args.prnet_source,args.transigen_source)):parser.error('explicit inputs, means and both source paths required')
    if args.import_cold_prnet:
        if any(p is None for p in (args.import_original,args.import_runtime,args.compatibility_audit)):
            parser.error('import requires original candidate, frozen runtime and reviewed compatibility audit')
        for run in plan():
            if run['adapter']=='prnet' and run['scenario']=='cold-drug':
                target=run_path(args.candidate,run)
                if args.resume and target.exists():verify_execution(target,run,frozen_plan(args))
                else:import_run(args,run)
    elif args.all:execute(args)
    else:accept(args)

if __name__=='__main__':main()
