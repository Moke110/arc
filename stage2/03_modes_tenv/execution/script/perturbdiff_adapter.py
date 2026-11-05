"""Native PerturbDiff with jointly trained Morgan input at its chemical-vector seam."""
import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import pickle
import random
import resource
import signal
import sys
import time
from types import SimpleNamespace
sys.dont_write_bytecode = True
for name in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(name,'2')
import numpy as np
import torch
from torch import nn
from omegaconf import OmegaConf
from artifacts import read_json, write_json, sha256
from predictor_inputs import PredictorInputs

REPO=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
UPSTREAM=REPO/'pert_model/PerturbDiff/source'
sys.path.insert(0,str(UPSTREAM))
from src.models.lightning.lightning_module import PlModel
from src.models.lightning.lightning_factories import get_optimizer, EMAWeightAveraging
from src.common.utils import get_cosine_schedule_with_warmup


def digest(module):
    result=hashlib.sha256()
    for name,value in module.state_dict().items():
        result.update(name.encode())
        result.update(value.detach().cpu().numpy().tobytes())
    return result.hexdigest()


def save(path,value):
    tmp=path.with_suffix(path.suffix+'.partial')
    torch.save(value,tmp)
    tmp.replace(path)


@contextmanager
def validation_random_stream(seed,device):
    state=np.random.get_state()
    try:
        np.random.seed(seed)
        with torch.random.fork_rng(devices=[device.index or 0] if device.type=='cuda' else []):
            torch.manual_seed(seed)
            yield
    finally:
        np.random.set_state(state)


class JointPerturbDiff(nn.Module):
    def __init__(self,inputs,config,directory):
        super().__init__()
        all_cases=[r for cases in inputs.cases.values() for r in cases]
        self.expression_scale=config['expression_scale']
        self.keys=sorted({(r['drug_key'],float(r['dose_um'])) for r in all_cases})
        self.lookup={key:i for i,key in enumerate(self.keys)}
        # Native CovEncoder builds its categorical dose vocabulary from metadata.
        # Temporary zero vectors establish the native vector shape; no learned
        # external embeddings enter training. The table is immediately replaced.
        native_vectors={key[0]:np.zeros(config['chemical_dimension'],dtype='float32') for key in self.keys}
        vector_path=directory/'native_vector_shapes.pkl'
        with vector_path.open('wb') as stream:
            pickle.dump(native_vectors,stream)
        model_cfg=OmegaConf.load(UPSTREAM/'configs/model/base.yaml')
        model_cfg.hidden_num=[inputs.n_genes,512]
        model_cfg.input_dim=inputs.n_genes
        model_cfg.dataset_dict={'tahoe_tenv':0}
        model_cfg.p_drop_control=0.0  # native scratch recipe, README COMMON_SCRATCH_DATA
        model_cfg.class_emb_gather_strategy='concat'
        model_cfg.class_emb_hidden_dimension=128
        cov_cfg=OmegaConf.create({'pert_encoding':'onehot','drug_encoding':'chemberta_cls',
            'replogle_gene_encoding':'onehot','hidden_dim':128,'output_dim':inputs.n_genes,
            'pert_dict':{repr([(key[0],key[1],'uM')]):i for i,key in enumerate(self.keys)},
            'num_pert':len(self.keys),'celltype_encoding':'onehot','num_celltype':1,'batch_encoding':None,
            'drug_embedding_path':str(vector_path)})
        # Native loader expects its zero-dose control entry in the metadata.
        cov_cfg.pert_dict[repr([('DMSO_TF',0.0,'uM')])]=len(self.keys)
        optimizer_cfg=OmegaConf.create({'micro_batch_size':config['batch_groups'],
            'blur':.05,'use_mse_loss':True,'MMD_loss_factor':1.,
            'optimizer':{'adam_w_mode':True,'lr':config['learning_rate'],'eps':1e-8,
                         'betas':[.9,.98],'weight_decay':config['weight_decay']}})
        self.native=PlModel(cov_cfg,model_cfg,None,optimizer_cfg,
                            OmegaConf.create({'log_every_n_steps':100}),['validation'])
        native_projection=self.native.cov_encoder.drug_encoder[1]
        self.adapter=nn.Linear(2048,config['chemical_dimension'])
        fingerprints=np.stack([np.array(inputs.morgan[inputs.drug_index[key[0]]],dtype='float32') for key in self.keys]
                              +[np.zeros(2048,dtype='float32')])
        # Preserve the native drug projection, dose encoder, gather and backbone.
        self.native.cov_encoder.drug_encoder=nn.Sequential(
            nn.Embedding.from_pretrained(torch.from_numpy(fingerprints),freeze=True),
            self.adapter,native_projection)
        vector_path.unlink()

    def batch(self,inputs,role,groups,device,rng=None,population_size=None):
        indices=[]
        for group in groups:
            selected=np.asarray(group,dtype=int)
            if population_size is not None:
                selected=rng.choice(selected,size=population_size,replace=len(selected)<population_size)
            indices.append(selected)
        index_array=np.stack(indices)
        baseline,_,_,*target=inputs.batch(role,index_array.reshape(-1),device)
        size=index_array.shape
        cov=torch.tensor([self.lookup[(inputs.cases[role][int(group[0])]['drug_key'],
                                     float(inputs.cases[role][int(group[0])]['dose_um']))] for group in index_array],device=device).unsqueeze(1)
        batch={'cont_emb':baseline.reshape(*size,inputs.n_genes)/self.expression_scale,'cov_pert':cov,
               'cov_celltype':torch.zeros(len(groups),1,dtype=torch.long,device=device),
               'cov_batch':torch.zeros(len(groups),1,dtype=torch.long,device=device),'ds_name':['tahoe_tenv']*len(groups)}
        if target:
            batch['pert_emb']=target[0].reshape(*size,inputs.n_genes)/self.expression_scale
        return batch,index_array

    def sample(self,batch,config,noise):
        condition={'batch_emb':self.native._encode_covariates(batch),'cont_emb':batch['cont_emb'],
                   'gene_emb':None,'ds_name':batch['ds_name']}
        result,_=self.native.diffusion.ddim_sample_loop(self.native.model,tuple(noise.shape),
            self_condition=condition,noise=noise,device=noise.device,clip_denoised=True,
            start_time=config['sampling_start_time'],eta=0.,guidance_strength=1.)
        return result


def groups_for(cases):
    groups={}
    for i,case in enumerate(cases):
        groups.setdefault(case['group_id'],[]).append(i)
    return list(groups.values())


def fit_predict(request,config):
    config=dict(config)
    config['batch_groups']=max(1,int(config['batch_groups']*float(os.environ.get('TENV_BATCH_SCALE','1'))))
    if (config['max_steps']<1 or (config['validation_interval'] is not None and config['validation_interval']<1) or config['validation_max_batches']<1
        or config['population_size']<1 or not 1<=config['sampling_start_time']<=1000
        or config['threads'] not in (1,2)):
        raise ValueError('invalid bounded native PerturbDiff recipe')
    torch.set_num_threads(config['threads'])
    seed=int(request['identity']['seed'])
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    device=torch.device(('cuda' if torch.cuda.is_available() else 'cpu') if config['device']=='auto' else config['device'])
    if device.type=='cuda':
        torch.cuda.manual_seed_all(seed); torch.cuda.reset_peak_memory_stats()
    directory=Path(request['output_directory']); directory.mkdir(parents=True,exist_ok=True)
    inputs=PredictorInputs(request)
    if config['validation_interval'] is None:
        config['validation_interval']=max(1,(len(groups_for(inputs.cases['train']))+config['batch_groups']-1)//config['batch_groups'])
    model=JointPerturbDiff(inputs,config,directory).to(device)
    initial_adapter=model.adapter.weight.detach().cpu().clone()
    initial_predictor=digest(model.native.model)
    decay=[]; no_decay=[]
    for name,p in model.named_parameters():
        if p.requires_grad:
            (no_decay if 'bias' in name or 'norm' in name else decay).append(p)
    optimizer=get_optimizer([{'params':decay,'weight_decay':config['weight_decay']},
                             {'params':no_decay,'weight_decay':0.}],model.native.optimizer_cfg.optimizer)
    scheduler=get_cosine_schedule_with_warmup(optimizer,num_warmup_steps=config['warmup_steps'],
                    num_training_steps=config['schedule_training_steps'],min_ratio=.1,plateau_ratio=.1)
    ema=EMAWeightAveraging(decay=config['ema_decay'],update_steps=config['ema_update_steps'])
    context=SimpleNamespace(global_step=0,current_epoch=0)
    ema.setup(context,model.native,'fit')
    rng=np.random.default_rng(seed)
    groups={role:groups_for(cases) for role,cases in inputs.cases.items()}
    start=time.monotonic(); best=float('inf'); best_step=None; history=[]; step_seconds=[]; stopped=False
    config_hash=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    def stop(signum,frame):
        nonlocal stopped
        stopped=True
    signal.signal(signal.SIGTERM,stop); signal.signal(signal.SIGINT,stop)
    def progress(phase,**extra):
        write_json(directory/'progress.json',{'phase':phase,'training_seconds':time.monotonic()-start,**extra})
    def snapshot(step):
        return {'identity':request['identity'],'config':config,'config_sha256':config_hash,
                'model_state':model.state_dict(),'optimizer_state':optimizer.state_dict(),
                'scheduler_state':scheduler.state_dict(),'step':step,'best_validation_loss':best,
                'training_seconds':time.monotonic()-start,'torch_rng':torch.get_rng_state(),'numpy_rng':rng.bit_generator.state}
    progress('training',step=0)
    stop_reason='max_steps'
    step=0
    while step<config['max_steps']:
        order=rng.permutation(len(groups['train']))
        for offset in range(0,len(order),config['batch_groups']):
            if stopped or time.monotonic()-start>=config['max_training_seconds']:
                stop_reason='signal' if stopped else 'training_budget'
                break
            if device.type=='cuda':
                torch.cuda.synchronize()
            step_start=time.monotonic()
            model.train()
            batch,_=model.batch(inputs,'train',[groups['train'][i] for i in order[offset:offset+config['batch_groups']]],
                                device,rng,config['population_size'])
            optimizer.zero_grad(set_to_none=True)
            losses=model.native._compute_loss(batch)
            loss=(losses['loss']*losses['weights']).mean()
            if not torch.isfinite(loss):
                raise ValueError('nonfinite native objective')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(),1.0)
            optimizer.step(); scheduler.step(); step+=1
            context.global_step=step
            ema.on_train_batch_end(context,model.native,None,batch,offset//config['batch_groups'])
            if device.type=='cuda':
                torch.cuda.synchronize()
            step_seconds.append(time.monotonic()-step_start)
            progress('training',step=step)
            if step%config['validation_interval']==0 or step==config['max_steps']:
                model.eval(); values=[]
                ema.on_validation_epoch_start(context,model.native)
                # Native validation computes the same randomly sampled diffusion
                # objective; seed its stream to compare checkpoints consistently.
                valid_rng=np.random.default_rng(seed+100000)
                with validation_random_stream(seed+100000,device):
                    with torch.no_grad():
                        for off in range(0,min(len(groups['validation']),config['batch_groups']*config['validation_max_batches']),config['batch_groups']):
                            vb,_=model.batch(inputs,'validation',groups['validation'][off:off+config['batch_groups']],
                                             device,valid_rng,config['population_size'])
                            vl=model.native._compute_loss(vb)
                            values.append(float((vl['loss']*vl['weights']).mean()))
                validation=float(np.mean(values))
                if not np.isfinite(validation):
                    raise ValueError('nonfinite native validation objective')
                history.append({'step':step,'training_loss':float(loss.detach()),'validation_loss':validation})
                if validation<best:
                    best=validation; best_step=step; save(directory/'best.pt',snapshot(step))
                ema.on_validation_epoch_end(context,model.native)
                progress('training',step=step,best_validation_loss=best)
            if step>=config['max_steps']:
                break
        ema.on_train_epoch_end(context,model.native)
        context.current_epoch+=1
        if stop_reason!='max_steps':
            break
    if stopped or stop_reason=='training_budget':
        save(directory/'interrupted.pt',snapshot(step))
        raise RuntimeError('PerturbDiff interrupted; partial output is not completion')
    if best_step is None:
        save(directory/'interrupted.pt',snapshot(step))
        raise RuntimeError('PerturbDiff budget ended before native validation checkpoint')
    training_seconds=time.monotonic()-start
    inference_start=time.monotonic()
    progress('inference')
    checkpoint=torch.load(directory/'best.pt',map_location=device,weights_only=False)
    if checkpoint['identity']!=request['identity'] or checkpoint['config_sha256']!=config_hash:
        raise ValueError('checkpoint identity mismatch')
    model.load_state_dict(checkpoint['model_state']); del checkpoint
    model.eval()
    probe,_=model.batch(inputs,'test',[groups['test'][0]],device)
    probe_noise=torch.zeros_like(probe['cont_emb'])
    with torch.no_grad():
        original=model.sample(probe,config,probe_noise)
        saved_control=probe['cont_emb']; probe['cont_emb']=torch.zeros_like(saved_control)
        baseline_sensitivity=float(torch.linalg.vector_norm(original-model.sample(probe,config,probe_noise)))
        probe['cont_emb']=saved_control
        dose_encoder=model.native.cov_encoder.dose_encoder
        original_dose=dose_encoder.weight.detach().clone()
        dose_encoder.weight.zero_()
        dose_sensitivity=float(torch.linalg.vector_norm(original-model.sample(probe,config,probe_noise)))
        dose_encoder.weight.copy_(original_dose)
        adapter_change=float(torch.linalg.vector_norm(model.adapter.weight.detach().cpu()-initial_adapter))
        predictor_changed=digest(model.native.model)!=initial_predictor
        model.load_state_dict(torch.load(directory/'best.pt',map_location=device,weights_only=False)['model_state'])
        reload_error=float((original-model.sample(probe,config,probe_noise)).abs().max())
    if not adapter_change>0 or not predictor_changed or not baseline_sensitivity>0 or not dose_sensitivity>0 or reload_error!=0:
        raise ValueError('native adapter/baseline/dose/reload gate failed')
    prediction=np.lib.format.open_memmap(directory/'prediction.npy.partial',mode='w+',dtype='float32',
                                          shape=(len(inputs.cases['test']),inputs.n_genes))
    # Batch populations of equal well count; never mix drug-dose conditions.
    buckets={}
    for group in groups['test']:
        buckets.setdefault(len(group),[]).append(group)
    generator=torch.Generator(device=device).manual_seed(seed+200000)
    with torch.no_grad():
        for population_groups in buckets.values():
            for offset in range(0,len(population_groups),config['batch_groups']):
                if stopped or time.monotonic()-inference_start>=config['inference_max_seconds']:
                    raise RuntimeError('PerturbDiff inference budget/interruption; partial is not completion')
                batch,index=model.batch(inputs,'test',population_groups[offset:offset+config['batch_groups']],device)
                noise=torch.randn(batch['cont_emb'].shape,device=device,generator=generator)
                delta=((model.sample(batch,config,noise)-batch['cont_emb'])*config['expression_scale']).cpu().numpy()
                if not np.isfinite(delta).all():
                    raise ValueError('nonfinite native full-axis prediction')
                prediction[index.reshape(-1)]=delta.reshape(-1,inputs.n_genes)
    prediction.flush(); del prediction
    (directory/'prediction.npy.partial').replace(directory/'prediction.npy')
    write_json(directory/'case_ids.json',[r['case_id'] for r in inputs.cases['test']])
    result={'identity':request['identity'],'author_objective':'native GaussianDiffusion MSE + energy MMD',
            'architecture':{'depth':8,'hidden_size':512,'diffusion_steps':1000},
            'best_step':best_step,'best_validation_loss':best,'history':history,'training_step_seconds':step_seconds,
            'stop_reason':stop_reason,
            'training_seconds':training_seconds,'inference_seconds':time.monotonic()-inference_start,
            'behavior':{'adapter_parameter_change_l2':adapter_change,'predictor_parameters_changed':predictor_changed,
                        'baseline_sensitivity_l2':baseline_sensitivity,'dose_sensitivity_l2':dose_sensitivity,
                        'reload_max_abs_error':reload_error},
            'trainable_parameters':sum(p.numel() for p in model.parameters() if p.requires_grad),
            'adapter_parameters':sum(p.numel() for p in model.adapter.parameters()),
            'blobs':{'best.pt':sha256(directory/'best.pt')},'config':config,'config_sha256':config_hash,
            'peak_rss_kbytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'peak_gpu_bytes':torch.cuda.max_memory_allocated() if device.type=='cuda' else 0,
            'peak_gpu_reserved_bytes':torch.cuda.max_memory_reserved() if device.type=='cuda' else 0,
            'environment':{'python':sys.version,'torch':torch.__version__,'numpy':np.__version__,'device':str(device)}}
    write_json(directory/'checkpoint.json',result)
    write_json(directory/'progress.json',{'phase':'complete','training_seconds':training_seconds,
                                        'inference_seconds':result['inference_seconds']})
    print(json.dumps({key:result[key] for key in ('best_step','best_validation_loss','behavior','training_seconds','inference_seconds')}),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request',required=True)
    parser.add_argument('--config',default=str(REPO/'stage2/03_modes_tenv/execution/models/perturbdiff/config.json'))
    args=parser.parse_args()
    fit_predict(read_json(args.request),read_json(args.config))
