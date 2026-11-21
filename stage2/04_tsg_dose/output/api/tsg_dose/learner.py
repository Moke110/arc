"""Validation-first native fit and separate post-choice test/replay lifecycle.

Reuse accepted02 native adapter, strict stopper and isolated RNG. No fit reads
outer/held-dose expression or emits test predictions before a choice receipt.
"""
import json
from functools import wraps
import os
import platform
import resource
import secrets
import time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from . import shared
from tenv5k.learning import EarlyStopping, fixed_rng
from tenv5k.attempt_identity import attempt_binding
from .inputs import Inputs, sha256
from .protocol import METHODS, SEEDS, score_prediction, summary
from .comparison import choose_method
from . import model as models
PROCESS_TOKEN=secrets.token_hex(32)


class ScientificFailure(ValueError):
    """Numerical or exact-replay failure that cannot be resumed or rescued."""


def scientific_phase(phase):
    def decorate(function):
        @wraps(function)
        def run(inputs,candidate,*args,**kwargs):
            try:return function(inputs,candidate,*args,**kwargs)
            except ScientificFailure as error:
                candidate=Path(candidate)
                terminal=json.loads((candidate/'terminal.json').read_text())
                evidence=dict(phase=phase,error=str(error),previous_terminal=terminal)
                if (candidate/'manifest.json').exists():evidence['previous_manifest']=json.loads((candidate/'manifest.json').read_text())
                write_json(candidate/'failure.json',evidence)
                write_json(candidate/'terminal.json',dict(terminal,status='failed',failure_phase=phase,error=str(error),previous_terminal=terminal))
                manifest(candidate,inputs)
                raise
        return run
    return decorate


def write_json(path,value):
    Path(path).write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')


def runtime_identity():
    """Owned API and exact reused mechanisms; no candidate-dependent hash scans."""
    here=Path(__file__).parent
    paths=list(here.glob('*.py'))
    paths += [shared.SHARED/'tenv5k'/name for name in ['model_transigen.py','native_import.py','learning.py','protocol.py','inputs.py','attempt_identity.py']]
    return {str(p.relative_to(here) if p.parent==here else Path('shared')/p.name):sha256(p) for p in sorted(paths)}


def environment(device):
    return dict(python=platform.python_version(),torch=torch.__version__,numpy=np.__version__,device=str(device),cuda=torch.version.cuda)


def serialize_prediction(value):
    with np.errstate(over='ignore',invalid='ignore'):value=np.asarray(value,dtype=np.float32)
    if not np.isfinite(value).all():raise ScientificFailure('non-finite serialized prediction')
    return value


def predict(native,model,data,rows,seed):
    with fixed_rng(seed,next(model.parameters()).device):return serialize_prediction(native.predict(model,data,rows))


def manifest(candidate,inputs):
    files={p.name:sha256(p) for p in sorted(Path(candidate).iterdir()) if p.is_file() and p.name!='manifest.json'}
    write_json(Path(candidate)/'manifest.json',dict(schema='tsg-dose-fit-v1',input_manifest_sha256=sha256(Path(inputs)/'manifest.json'),execution_binding=json.loads((Path(candidate)/'config.json').read_text())['execution_binding'],outputs=files))


def verify_candidate(candidate):
    candidate=Path(candidate);registered=json.loads((candidate/'manifest.json').read_text())
    if registered.get('schema')!='tsg-dose-fit-v1':raise ValueError('fit manifest schema mismatch')
    actual={p.name for p in candidate.iterdir() if p.is_file() and p.name!='manifest.json'}
    if actual!=set(registered['outputs']):raise ValueError('fit manifest coverage mismatch')
    for name,digest in registered['outputs'].items():
        if Path(name).name!=name or sha256(candidate/name)!=digest:raise ValueError('fit digest mismatch '+name)
    return registered


def smoke_subset(data,rows,limit):
    if limit is None:return rows
    buckets=[list(rows[data.doses(rows)==dose]) for dose in sorted(set(data.doses(rows)))]
    chosen=[]
    for level in range(max(map(len,buckets))):
        for bucket in buckets:
            if level<len(bucket):chosen.append(bucket[level])
            if len(chosen)==limit:return np.array(sorted(chosen),dtype=int)
    return np.array(sorted(chosen),dtype=int)


def train(inputs,candidate,source,*,method='no-dose',scenario='cold-drug',fold=0,seed=17,device='cpu',smoke_rows=None,max_epochs=1000,registration=None):
    candidate=Path(candidate)
    if candidate.exists() and any(candidate.iterdir()):raise ValueError('candidate must be empty')
    candidate.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    try:
        if method not in METHODS or scenario not in ['cold-drug','held-dose','final'] or fold not in (range(3) if scenario!='final' else [0]) or seed not in SEEDS:raise ValueError('unregistered fit identity')
        if smoke_rows is None and max_epochs!=1000:raise ValueError('formal fits require1000 epoch budget')
        if smoke_rows is not None and smoke_rows<2:raise ValueError('smoke requires at least2 rows')
        if registration is None:
            if smoke_rows is None:raise ValueError('formal fit requires preregistered execution attempt')
            registration=dict(attempt_nonce=secrets.token_hex(32),run=dict(method=method,scenario=scenario,fold=fold,seed=seed),device=str(device),origin='bounded-smoke')
        if registration.get('run')!=dict(method=method,scenario=scenario,fold=fold,seed=seed) or registration.get('device')!=str(device):raise ValueError('attempt registration identity mismatch')
        binding=attempt_binding(registration)
        write_json(candidate/'registration.json',registration)
        data=Inputs(inputs)
        train_rows=data.rows(scenario,fold,'train');valid_rows=data.rows(scenario,fold,'validation')
        if smoke_rows is not None:train_rows=smoke_subset(data,train_rows,smoke_rows);valid_rows=smoke_subset(data,valid_rows,smoke_rows)
        if not len(train_rows) or not len(valid_rows):raise ValueError('empty fitting partition')
        transform=models.fit_transform(method,data.doses(train_rows))
        native=models.adapter(source,method,transform)
        os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8');torch.use_deterministic_algorithms(True)
        backend=models.configure_backend(device)
        config=dict(models.configuration(method,transform),scenario=scenario,fold=fold,seed=seed,device=str(device),smoke_rows=smoke_rows,max_epochs=max_epochs,patience=50,validation_seed=22617,inference_seed=seed+100000,native_source_sha256=models.SOURCE_DIGEST,runtime_identity=runtime_identity(),backend=backend,train_order=train_rows.tolist(),validation_order=valid_rows.tolist(),execution_binding=binding)
        write_json(candidate/'config.json',config)
        model=native.build(seed,device);optimizer=native.optimizer(model)
        stopper=EarlyStopping(max_epochs=max_epochs);shuffle=np.random.default_rng(seed);history=[]
        initial={k:v.detach().cpu().clone() for k,v in model.state_dict().items() if k.endswith('weight')}
        for epoch in range(1,max_epochs+1):
            loss=native.train_epoch(model,data,train_rows,optimizer,shuffle)
            validation=predict(native,model,data,valid_rows,config['validation_seed'])
            mse=float(np.mean(np.square(validation-np.asarray(data.expression[valid_rows],dtype=np.float32))))
            improved=stopper.observe(mse,epoch)
            history.append(dict(epoch=epoch,train_loss=loss,validation_mse=mse,best_mse=stopper.best,stale=stopper.stale,lr=optimizer.param_groups[0]['lr']))
            if improved:torch.save(dict(state_dict=model.state_dict(),config=config,best_epoch=epoch,validation_mse=mse,input_manifest_sha256=sha256(Path(inputs)/'manifest.json')),candidate/'checkpoint.pt')
            if stopper.stopped:break
        if not any(not torch.equal(v,model.state_dict()[k].detach().cpu()) for k,v in initial.items()):raise ValueError('training did not update weights')
        saved=torch.load(candidate/'checkpoint.pt',map_location=device,weights_only=False);model.load_state_dict(saved['state_dict'])
        validation=predict(native,model,data,valid_rows,config['validation_seed'])
        np.save(candidate/'validation_prediction.npy',validation)
        wells,groups=score_prediction(data.obs.iloc[valid_rows],data.expression[valid_rows],validation,data.control[valid_rows],prediction_ids=data.obs.index[valid_rows])
        wells.to_csv(candidate/'validation_per_well.csv',index=False);groups.to_csv(candidate/'validation_per_group.csv',index=False)
        data.obs.iloc[valid_rows][['profile_id','drug_parent_inchikey','dose_um']].to_csv(candidate/'validation_identity.csv',index=False)
        write_json(candidate/'validation_metrics.json',summary(groups))
        pd.DataFrame(history).to_csv(candidate/'learning_curve.csv',index=False)
        write_json(candidate/'environment.json',environment(device))
        terminal=dict(training_process_token=PROCESS_TOKEN,training_pid=os.getpid(),status='validation-trained-awaiting-replay',smoke=smoke_rows is not None,best_epoch=stopper.best_epoch,epochs=len(history),validation_mse=stopper.best,elapsed_seconds=time.monotonic()-started,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,weights_updated=True,parameter_count=sum(p.numel() for p in model.parameters()))
        write_json(candidate/'terminal.json',terminal)
        return terminal
    except Exception as error:
        write_json(candidate/'terminal.json',dict(status='failed',error=str(error)));raise


def load(inputs,candidate,source):
    data=Inputs(inputs);candidate=Path(candidate)
    config=json.loads((candidate/'config.json').read_text())
    if config['execution_binding']!=attempt_binding(json.loads((candidate/'registration.json').read_text())):raise ValueError('checkpoint attempt registration mismatch')
    saved=torch.load(candidate/'checkpoint.pt',map_location='cpu',weights_only=False)
    if saved['config']!=config or saved['input_manifest_sha256']!=sha256(Path(inputs)/'manifest.json'):raise ValueError('checkpoint configuration/input mismatch')
    if config['runtime_identity']!=runtime_identity() or config['native_source_sha256']!=models.SOURCE_DIGEST:raise ValueError('checkpoint runtime/source mismatch')
    if json.loads((candidate/'environment.json').read_text())!=environment(config['device']):raise ValueError('replay environment mismatch')
    expected=models.fit_transform(config['method'],data.doses(np.asarray(config['train_order'],dtype=int)))
    if config['transform']!=expected or any(config[k]!=v for k,v in models.configuration(config['method'],expected).items()):raise ValueError('method/transform identity mismatch')
    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8');torch.use_deterministic_algorithms(True)
    if models.configure_backend(config['device'])!=config['backend']:raise ValueError('backend mismatch')
    if any(v.is_floating_point() and v.dtype!=torch.float32 for v in saved['state_dict'].values()):raise ValueError('checkpoint precision mismatch')
    native=models.adapter(source,config['method'],expected);model=native.build(config['seed'],config['device']);model.load_state_dict(saved['state_dict'])
    return data,native,model,config


@scientific_phase('replay')
def accept_replay(inputs,candidate,source):
    """Fresh process only: reproduce validation and, if present, test predictions."""
    candidate=Path(candidate)
    terminal=json.loads((candidate/'terminal.json').read_text())
    if terminal.get('training_pid')==os.getpid() or terminal.get('training_process_token')==PROCESS_TOKEN:raise ValueError('replay must run in a fresh process')
    if terminal['status'] not in ['validation-trained-awaiting-replay','validation-complete','test-awaiting-replay','complete']:raise ValueError('cannot replay failed candidate')
    if (candidate/'manifest.json').exists():verify_candidate(candidate)
    data,native,model,config=load(inputs,candidate,source)
    for partition,seed in [('validation',config['validation_seed']),('test',config['inference_seed'])]:
        path=candidate/(partition+'_prediction.npy')
        if not path.exists():continue
        rows=np.asarray(config['validation_order']) if partition=='validation' else data.rows(config['scenario'],config['fold'],'test')
        if partition=='test' and config['smoke_rows'] is not None:rows=smoke_subset(data,rows,config['smoke_rows'])
        replayed=predict(native,model,data,rows,seed);stored=np.load(path)
        if stored.dtype!=np.float32 or not np.array_equal(stored,replayed):raise ScientificFailure('exact fresh-process '+partition+' replay differs')
    phase='test' if (candidate/'test_prediction.npy').exists() else 'validation'
    terminal[phase+'_replay']=dict(pid=os.getpid(),process_token=PROCESS_TOKEN,attempt_nonce=config['execution_binding']['attempt_nonce'])
    terminal.update(replay_pid=os.getpid(),status='complete' if (candidate/'test_prediction.npy').exists() else 'validation-complete',replay='exact fresh-process native sampled prediction')
    write_json(candidate/'terminal.json',terminal)
    if not (candidate/'validation_manifest.json').exists():
        write_json(candidate/'validation_manifest.json',dict(schema='tsg-dose-validation-v1',input_manifest_sha256=sha256(Path(inputs)/'manifest.json'),parameter_count=terminal['parameter_count'],outputs={p.name:sha256(p) for p in sorted(candidate.iterdir()) if p.is_file() and p.name not in ['manifest.json','terminal.json','validation_manifest.json']}))
    manifest(candidate,inputs);return terminal


@scientific_phase('evaluation')
def evaluate(inputs,candidate,source,selection=None):
    """Access outer labels only after all21 validation fits and a locked choice."""
    candidate=Path(candidate);verify_candidate(candidate)
    config=json.loads((candidate/'config.json').read_text())
    terminal=json.loads((candidate/'terminal.json').read_text())
    if terminal['status']!='validation-complete':raise ValueError('evaluation requires replayed validation checkpoint')
    if config['scenario']=='final':raise ValueError('final internal selection has no independent test')
    receipt=None
    if config['smoke_rows'] is None:
        if selection is None:raise ValueError('outer evaluation requires locked selection receipt')
        receipt=json.loads(Path(selection).read_text())
        if receipt.get('scenario')!=config['scenario'] or receipt.get('fold')!=config['fold'] or receipt.get('input_manifest_sha256')!=sha256(Path(inputs)/'manifest.json'):raise ValueError('selection receipt input/partition mismatch')
        identities=set();scores=[];current_origin=False
        if len(receipt.get('validation_runs',{}))!=21:raise ValueError('selection requires exactly21 validation origins')
        for path,digest in receipt.get('validation_runs',{}).items():
            path=Path(path)
            if not path.is_absolute():path=Path(selection).parent/path
            if sha256(path/'validation_manifest.json')!=digest:raise ValueError('selection validation manifest changed')
            frozen=json.loads((path/'validation_manifest.json').read_text())
            if frozen.get('input_manifest_sha256')!=sha256(Path(inputs)/'manifest.json'):raise ValueError('validation input freeze mismatch')
            for name,expected in frozen['outputs'].items():
                if Path(name).name!=name or sha256(path/name)!=expected:raise ValueError('validation artifact changed')
            run=json.loads((path/'config.json').read_text());status=json.loads((path/'terminal.json').read_text())
            if run['scenario']!=config['scenario'] or run['fold']!=config['fold'] or run['smoke_rows'] is not None or status['status'] not in ['validation-complete','test-awaiting-replay','complete']:raise ValueError('invalid validation origin')
            identities.add((run['method'],run['seed']))
            parameter_count=frozen.get('parameter_count')
            if parameter_count is None or parameter_count!=status.get('parameter_count'):raise ValueError('validation parameter count changed')
            scores.append(dict(method=run['method'],seed=run['seed'],parameter_count=parameter_count,pcc=json.loads((path/'validation_metrics.json').read_text())['pcc']))
            if path.resolve()==candidate.resolve() and digest==sha256(candidate/'validation_manifest.json'):current_origin=True
        if identities!={(m,s) for m in METHODS for s in SEEDS}:raise ValueError('selection requires21 replayed formal validation fits')
        if not current_origin:raise ValueError('evaluated candidate is absent from locked validation origins')
        chosen=choose_method(pd.DataFrame(scores))
        expected='no-dose' if config['scenario']=='held-dose' else chosen
        if receipt.get('method')!=expected:raise ValueError('selection does not match validation winner')
    data,native,model,config=load(inputs,candidate,source)
    rows=data.rows(config['scenario'],config['fold'],'test')
    if config['smoke_rows'] is not None:rows=smoke_subset(data,rows,config['smoke_rows'])
    prediction=predict(native,model,data,rows,config['inference_seed']);np.save(candidate/'test_prediction.npy',prediction)
    data.obs.iloc[rows][['profile_id','drug_parent_inchikey','dose_um']].to_csv(candidate/'test_identity.csv',index=False)
    wells,groups=score_prediction(data.obs.iloc[rows],data.expression[rows],prediction,data.control[rows],prediction_ids=data.obs.index[rows])
    wells.to_csv(candidate/'test_per_well.csv',index=False);groups.to_csv(candidate/'test_per_group.csv',index=False)
    if selection is not None:write_json(candidate/'selection_receipt.json',dict(receipt=receipt,sha256=sha256(selection)))
    terminal.update(status='test-awaiting-replay');write_json(candidate/'terminal.json',terminal)
    # Existing validation manifest remains sealed in selection. Preserve it and
    # replace the candidate manifest only after exact test replay in a child.
    (candidate/'manifest.json').unlink()
    return terminal
