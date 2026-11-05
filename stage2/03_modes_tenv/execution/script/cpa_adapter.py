"""Native CPA module/training plan, streamed frozen inputs and joint Morgan map."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import random
import resource
import signal
import sys
import time

for name in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(name,'2')
os.environ.setdefault('MPLCONFIGDIR','/tmp/arc-cpa-matplotlib')
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/arc-cpa-numba')
import numpy as np
import torch
from torch import nn
from artifacts import read_json,write_json,sha256
from predictor_inputs import PredictorInputs,minibatches


def repo_root():
    for candidate in Path(__file__).resolve().parents:
        if (candidate/'.git').exists(): return candidate
    raise RuntimeError('repository root not found')

ROOT=repo_root()
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'pert_model/CPA/source'))
# Source checkout is not installed as a distribution; match the setup run's fallback.
_original_version=importlib.metadata.version

def _checkout_version(name):
    try: return _original_version(name)
    except importlib.metadata.PackageNotFoundError:
        if name=='cpa-tools': return '0.0.0+source'
        raise

importlib.metadata.version=_checkout_version
try:
    from cpa._module import CPAModule
    from cpa._task import CPATrainingPlan
    from cpa._utils import CPA_REGISTRY_KEYS
finally:
    importlib.metadata.version=_original_version


class StreamingPlan(CPATrainingPlan):
    """Use native training_step with an explicit optimizer/backward host.

    This avoids materializing a fold AnnData or introducing a Lightning loader.
    Losses, adversarial penalties and alternating steps remain upstream calls.
    """
    @property
    def current_epoch(self): return getattr(self,'stream_epoch',0)
    @property
    def global_step(self): return getattr(self,'stream_step',0)
    def optimizers(self): return self.stream_optimizers
    def manual_backward(self,loss,*args,**kwargs): loss.backward(*args,**kwargs)


def vector(module):
    return torch.cat([p.detach().cpu().reshape(-1) for p in module.parameters() if p.requires_grad])


def save_torch(path,payload):
    temp=path.with_suffix(path.suffix+'.partial')
    torch.save(payload,temp);temp.replace(path)


def make_plan(inputs,n_genes,config,seed,device):
    fingerprints=torch.zeros(len(inputs.drug_index)+1,2048)
    for drug,i in inputs.drug_index.items():
        fingerprints[i+1]=torch.from_numpy(np.array(inputs.morgan[i],dtype='float32'))
    embedding=nn.Embedding.from_pretrained(fingerprints,freeze=True,padding_idx=0)
    module=CPAModule(n_genes=n_genes,n_perts=len(fingerprints),covars_encoder={},
        drug_embeddings=embedding,n_latent=config['n_latent'],recon_loss='gauss',doser_type='linear',
        n_hidden_encoder=config['n_hidden_encoder'],n_hidden_decoder=config['n_hidden_decoder'],
        n_layers_encoder=config['n_layers_encoder'],n_layers_decoder=config['n_layers_decoder'],
        dropout_rate_encoder=config['dropout_rate_encoder'],dropout_rate_decoder=config['dropout_rate_decoder'],
        variational=False,seed=seed)
    plan=StreamingPlan(module,covars_to_ncovars={},n_adv_perts=len(fingerprints),
        lr=config['learning_rate'],wd=config['weight_decay'],adv_lr=config['adv_learning_rate'],
        adv_wd=config['adv_weight_decay'],adv_steps=config['adv_steps'],reg_adv=config['reg_adv'],
        pen_adv=config['pen_adv'],n_hidden_adv=config['n_hidden_adv'],n_layers_adv=config['n_layers_adv'],
        step_size_lr=config['step_size_lr'],mixup_alpha=0.)
    plan=plan.to(device)
    optimizers,schedulers=plan.configure_optimizers()
    # Upstream configure_optimizers omits the external chem projection. The agreed
    # adapter must train jointly, so attach it to the native autoencoder optimizer.
    optimizers[0].add_param_group({'params':module.pert_network.pert_transformation.parameters()})
    plan.stream_optimizers=optimizers
    return plan,schedulers


def native_batch(inputs,role,indices,device):
    """Author reconstruction on observed expression plus allowed control rows."""
    count=len(inputs.cases[role])
    expressed=[];pert=[];dose=[];category=[]
    group_categories={r['group_id']:i for i,r in enumerate(inputs.cases[role])}
    for index in indices:
        index=int(index)
        if index<count:
            case=inputs.cases[role][index]
            control=int(case['control_index'])
            expressed.append(np.array(inputs.delta[role][index],dtype='float32')+inputs.controls[control])
            pert.append(inputs.drug_index[case['drug_key']]+1)
            dose.append(float(case['dose_um']))
            category.append(group_categories[case['group_id']])
        else:
            control=index-count
            expressed.append(np.array(inputs.controls[control],dtype='float32'))
            pert.append(0);dose.append(0.);category.append(index)
    x=torch.tensor(np.array(expressed),device=device)
    if not torch.isfinite(x).all(): raise ValueError('nonfinite native CPA training inputs')
    return {CPA_REGISTRY_KEYS.X_KEY:x,
            CPA_REGISTRY_KEYS.PERTURBATION_KEY:torch.tensor(pert,device=device),
            CPA_REGISTRY_KEYS.PERTURBATIONS:torch.tensor(pert,device=device)[:,None],
            CPA_REGISTRY_KEYS.PERTURBATIONS_DOSAGES:torch.tensor(dose,device=device)[:,None],
            CPA_REGISTRY_KEYS.CATEGORY_KEY:torch.tensor(category,device=device)}


def predict(module,inputs,indices,device,zero_baseline=False):
    cases=[inputs.cases['test'][int(i)] for i in indices]
    baseline=torch.from_numpy(np.array(inputs.controls[[int(r['control_index']) for r in cases]],dtype='float32')).to(device)
    perts=torch.tensor([inputs.drug_index[r['drug_key']]+1 for r in cases],device=device)[:,None]
    doses=torch.tensor([float(r['dose_um']) for r in cases],device=device)[:,None]
    inference=module.inference(x=torch.zeros_like(baseline) if zero_baseline else baseline,
             perts={'true':perts},perts_doses={'true':doses},covars_dict={})
    expression=module.generative(z=inference['z'],library=inference['library'])['px'].loc
    return expression,baseline


def fit_predict(request,config):
    config=dict(config)
    scale=float(os.environ.get('TENV_BATCH_SCALE','1'))
    if not 0<scale<=1: raise ValueError('invalid retry batch scale')
    config['batch_size']=max(2,int(config['batch_size']*scale))
    if config['batch_size']<2 or config['max_epochs']<1 or config['threads'] not in (1,2):
        raise ValueError('invalid bounded CPA recipe')
    CPA_REGISTRY_KEYS.PERTURBATION_KEY='perturbation'
    torch.set_num_threads(config['threads'])
    seed=int(request['identity']['seed'])
    random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
    device=torch.device(('cuda' if torch.cuda.is_available() else 'cpu') if config['device']=='auto' else config['device'])
    if device.type=='cuda': torch.cuda.reset_peak_memory_stats()
    inputs=PredictorInputs(request)
    output=Path(request['output_directory']);output.mkdir(parents=True,exist_ok=True)
    plan,schedulers=make_plan(inputs,request['n_genes'],config,seed,device)
    initial_adapter=vector(plan.module.pert_network.pert_transformation)
    initial_predictor=vector(plan.module)
    initial_adversary=vector(plan.perturbation_classifier)
    best=float('inf');best_epoch=None;missed=0;history=[]
    rng=np.random.default_rng(seed)
    start=time.monotonic();stopped=False;stop_reason='max_epochs'
    config_digest=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()

    def stop(signum,frame):
        nonlocal stopped
        stopped=True
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    def snapshot(epoch):
        return {'identity':request['identity'],'config_sha256':config_digest,'plan_state':plan.state_dict(),
                'optimizers':[o.state_dict() for o in plan.stream_optimizers],
                'schedulers':[s.state_dict() for s in schedulers],'epoch':epoch,
                'training_seconds':time.monotonic()-start,'config':config}
    def progress(phase,**extra):
        write_json(output/'progress.json',{'phase':phase,'training_seconds':time.monotonic()-start,**extra})
    train_count=len(inputs.cases['train'])+len(inputs.controls)
    for epoch in range(config['max_epochs']):
        plan.train();plan.stream_epoch=epoch
        total=0.;observations=0
        for batch_idx,indices in enumerate(minibatches(train_count,config['batch_size'],rng,training=True)):
            if stopped or time.monotonic()-start>=config['max_training_seconds']:
                stopped=True;break
            batch=native_batch(inputs,'train',indices,device)
            result=plan.training_step(batch,batch_idx)
            if not np.isfinite(result['recon_loss']): raise ValueError('nonfinite native CPA loss')
            total+=result['recon_loss']*len(indices);observations+=len(indices)
            plan.stream_step=getattr(plan,'stream_step',0)+1
            progress('training',epoch=epoch+1)
        if stopped:
            save_torch(output/'interrupted.pt',snapshot(epoch));progress('stopped',stop_reason='signal_or_training_budget')
            raise RuntimeError('CPA interrupted; saved state, no completed prediction')
        plan.eval();validation=0.;valid_count=0
        with torch.no_grad():
            for indices in minibatches(len(inputs.cases['validation']),config['batch_size']):
                batch=native_batch(inputs,'validation',indices,device)
                batch,lam=plan.module.mixup_data(batch,alpha=0.)
                inf,gen=plan.module.forward(batch,compute_loss=False,inference_kwargs={'mixup_lambda':lam})
                recon,_=plan.module.loss(batch,inf,gen)
                validation+=float(recon)*len(indices);valid_count+=len(indices)
        validation/=valid_count
        if stopped or time.monotonic()-start>=config['max_training_seconds']:
            save_torch(output/'interrupted.pt',snapshot(epoch));progress('stopped',stop_reason='signal_or_training_budget')
            raise RuntimeError('CPA interrupted during validation; no completed prediction')
        if not np.isfinite(validation): raise ValueError('nonfinite native CPA validation objective')
        history.append({'epoch':epoch+1,'training_loss':total/observations,'validation_loss':validation})
        if validation<best:
            best=validation;best_epoch=epoch+1;missed=0
            save_torch(output/'best.pt',snapshot(epoch+1))
        else: missed+=1
        # Match the native training_epoch_end StepLR cadence.
        if epoch>1 and epoch%config['step_size_lr']==0:
            for scheduler in schedulers: scheduler.step()
        progress('training',epoch=epoch+1,best_epoch=best_epoch,best_validation_loss=best)
        if missed>=config['early_stop_patience']:
            stop_reason='author_patience';break
    training_seconds=time.monotonic()-start
    if stopped:
        save_torch(output/'interrupted.pt',snapshot(epoch));progress('stopped',stop_reason='signal')
        raise RuntimeError('CPA interrupted during validation')
    inf_start=time.monotonic()
    write_json(output/'progress.json',{'phase':'inference','training_seconds':training_seconds})
    saved=torch.load(output/'best.pt',map_location=device)
    if saved['identity']!=request['identity'] or saved['config_sha256']!=config_digest:
        raise ValueError('stale CPA selected checkpoint')
    plan.load_state_dict(saved['plan_state']);plan.eval()
    probe=np.arange(min(2,len(inputs.cases['test'])))
    with torch.no_grad():
        original,_=predict(plan.module,inputs,probe,device)
        changed,_=predict(plan.module,inputs,probe,device,zero_baseline=True)
        behavior={'adapter_parameter_change_l2':float(torch.linalg.vector_norm(vector(plan.module.pert_network.pert_transformation)-initial_adapter)),
                  'predictor_parameter_change_l2':float(torch.linalg.vector_norm(vector(plan.module)-initial_predictor)),
                  'adversary_parameter_change_l2':float(torch.linalg.vector_norm(vector(plan.perturbation_classifier)-initial_adversary)),
                  'baseline_sensitivity_l2':float(torch.linalg.vector_norm(original-changed))}
        reloaded,_=make_plan(inputs,request['n_genes'],config,seed,device)
        reloaded.load_state_dict(torch.load(output/'best.pt',map_location=device)['plan_state']);reloaded.eval()
        replayed,_=predict(reloaded.module,inputs,probe,device)
        behavior['reload_max_abs_error']=float((original-replayed).abs().max())
    if any(behavior[k]<=0 for k in ['adapter_parameter_change_l2','predictor_parameter_change_l2','adversary_parameter_change_l2','baseline_sensitivity_l2']) or behavior['reload_max_abs_error']!=0.:
        raise ValueError('native CPA adapter/baseline/adversary/reload behavioral gate failed')
    del plan,initial_adapter,initial_predictor,initial_adversary
    prediction=np.lib.format.open_memmap(output/'prediction.npy.partial',mode='w+',dtype='float32',shape=(len(inputs.cases['test']),request['n_genes']))
    with torch.no_grad():
        for indices in minibatches(len(inputs.cases['test']),config['batch_size']):
            if stopped or time.monotonic()-inf_start>=config['inference_max_seconds']:
                raise RuntimeError('CPA incomplete inference; partial output not completion')
            expression,baseline=predict(reloaded.module,inputs,indices,device)
            delta=(expression-baseline).cpu().numpy()
            if not np.isfinite(delta).all(): raise ValueError('nonfinite CPA full-axis prediction')
            prediction[indices]=delta
    prediction.flush();del prediction
    (output/'prediction.npy.partial').replace(output/'prediction.npy')
    write_json(output/'case_ids.json',[r['case_id'] for r in inputs.cases['test']])
    inference_seconds=time.monotonic()-inf_start
    result={'identity':request['identity'],'author_objective':'CPA native Gaussian reconstruction NLL + native adversarial disentanglement',
            'best_epoch':best_epoch,'best_validation_loss':best,'stop_reason':stop_reason,'history':history,
            'training_seconds':training_seconds,'inference_seconds':inference_seconds,
            'selection':'minimum native Gaussian reconstruction NLL on validation',
            'trainable_parameters':sum(p.numel() for p in reloaded.parameters() if p.requires_grad),
            'adapter_parameters':sum(p.numel() for p in reloaded.module.pert_network.pert_transformation.parameters()),
            'behavior':behavior,'config':config,'config_sha256':config_digest,'blobs':{'best.pt':sha256(output/'best.pt')},
            'source':'CPA CPAModule + CPATrainingPlan.training_step; unchanged upstream',
            'recipe_deviations':['fixed full gene axis','shared Morgan external embeddings replace native extraction','native linear dose option replaces drug-specific learned dose scale','add upstream-omitted chem projection to AE optimizer','streamed native training_step replaces Lightning/AnnData materialization','no varying cell-type covariate because tenv has one declared context'],
            'environment':{'python':sys.version,'torch':torch.__version__,'numpy':np.__version__,'scvi_tools':_original_version('scvi-tools'),'device':str(device)},
            'peak_rss_kbytes':int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
            'peak_gpu_bytes':int(torch.cuda.max_memory_allocated()) if device.type=='cuda' else 0,
            'peak_gpu_reserved_bytes':int(torch.cuda.max_memory_reserved()) if device.type=='cuda' else 0}
    write_json(output/'checkpoint.json',result)
    write_json(output/'progress.json',{'phase':'complete','training_seconds':training_seconds,'inference_seconds':inference_seconds})
    print(json.dumps({'best_epoch':best_epoch,'best_validation_loss':best,'behavior':behavior}),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request',required=True)
    parser.add_argument('--config',default=str(Path(__file__).resolve().parent.parent/'models/cpa/config.json'))
    args=parser.parse_args()
    fit_predict(read_json(args.request),read_json(args.config))
