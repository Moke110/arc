"""Native ComPert with joint Morgan chemistry and native adversarial objectives."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import resource
import signal
import sys
import time

sys.dont_write_bytecode = True
for _name in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_name,'2')
import numpy as np
import torch
from torch import nn
from artifacts import read_json, write_json, sha256
from predictor_inputs import PredictorInputs, minibatches


def repo_root():
    for candidate in Path(__file__).resolve().parents:
        if (candidate/'.git').exists():
            return candidate
    raise RuntimeError('repository root not found')


sys.path.insert(0,str(repo_root()/'pert_model/chemCPA/source'))
from chemCPA.model import ComPert


class MorganRepresentation(nn.Module):
    """The native embedding-table interface, evaluated through one joint Linear."""
    def __init__(self, fingerprints, dimension):
        super().__init__()
        self.embedding_dim = dimension
        self.register_buffer('fingerprints',torch.tensor(np.array(fingerprints),dtype=torch.float32))
        self.linear = nn.Linear(2048,dimension)

    @property
    def weight(self):
        # Native ComPert gathers from .weight; no independent trainable drug rows.
        return self.linear(self.fingerprints)


def make_model(inputs, config, device, seed):
    table = MorganRepresentation(inputs.morgan,config['chemical_dimension'])
    model = ComPert(num_genes=inputs.n_genes,num_drugs=len(inputs.morgan),num_covariates=[0],
                    device=str(device),seed=seed,patience=config['early_stop_patience'],
                    doser_type='amortized',decoder_activation='ReLU',hparams=config['hparams'],
                    drug_embeddings=table,use_drugs_idx=True,enable_cpa_mode=False)
    return model.to(device)


def vector(modules):
    return torch.cat([p.detach().cpu().reshape(-1) for module in modules for p in module.parameters()])


def save_checkpoint(path, state):
    temporary = path.with_suffix(path.suffix+'.partial')
    torch.save(state,temporary)
    temporary.replace(path)


def context(inputs, role, indices, device):
    batch = inputs.batch(role,indices,device)
    drug_ids = torch.tensor([inputs.drug_index[inputs.cases[role][int(i)]['drug_key']] for i in indices],
                            dtype=torch.long,device=device)
    return batch,drug_ids


def forward(model, genes, drug_ids, doses):
    output,_,latents = model.predict(genes=genes,drugs_idx=drug_ids,dosages=doses,
                                    covariates=[],return_latent_basal=True)
    mean,variance = output.chunk(2,dim=1)
    return mean,variance,latents[0]


def losses(model, mean, variance, latent, expression, drug_ids):
    reconstruction = model.loss_autoencoder(mean,expression,variance)
    predictions = model.adversary_drugs(latent)
    target = torch.zeros(len(drug_ids),model.num_drugs,device=expression.device)
    target[torch.arange(len(drug_ids),device=expression.device),drug_ids] = 1.
    adversary = model.loss_adversary_drugs(predictions,target)
    autoencoder = reconstruction-model.hparams['reg_adversary']*adversary
    return reconstruction,adversary,autoencoder,predictions


def fit_predict(request, config):
    config = dict(config)
    batch_scale = float(os.environ.get('TENV_BATCH_SCALE','1'))
    if not 0 < batch_scale <= 1:
        raise ValueError('invalid scheduler smaller-batch scale')
    config['configured_batch_size'] = config['batch_size']
    config['batch_size'] = max(2,int(config['batch_size']*batch_scale))
    config['batch_scale'] = batch_scale
    if config['batch_size'] < 2 or config['max_epochs'] < 1 or config['threads'] not in (1,2):
        raise ValueError('invalid bounded chemCPA recipe')
    seed = int(request['identity']['seed'])
    torch.set_num_threads(config['threads'])
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.cuda.reset_peak_memory_stats()
    device = torch.device(('cuda' if torch.cuda.is_available() else 'cpu') if config['device']=='auto' else config['device'])
    directory = Path(request['output_directory'])
    directory.mkdir(parents=True,exist_ok=True)
    inputs = PredictorInputs(request)
    model = make_model(inputs,config,device,seed)
    ae_modules = [model.encoder,model.decoder,model.drug_embedding_encoder,model.drug_embeddings.linear]
    ae_parameters = [p for module in ae_modules for p in module.parameters()]
    initial = {'adapter':vector([model.drug_embeddings.linear]),
               'predictor':vector([model.encoder,model.decoder,model.drug_embedding_encoder]),
               'adversary':vector([model.adversary_drugs]),'doser':vector([model.dosers])}
    hp = model.hparams
    optimizers = [torch.optim.Adam(ae_parameters,lr=hp['autoencoder_lr'],weight_decay=hp['autoencoder_wd']),
                  torch.optim.Adam(model.adversary_drugs.parameters(),lr=hp['adversary_lr'],weight_decay=hp['adversary_wd']),
                  torch.optim.Adam(model.dosers.parameters(),lr=hp['dosers_lr'],weight_decay=hp['dosers_wd'])]
    schedulers = [torch.optim.lr_scheduler.StepLR(opt,step_size=hp['step_size_lr'],gamma=.9) for opt in optimizers]
    config_digest = hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    rng = np.random.default_rng(seed)
    history = []
    best = float('inf')
    best_epoch = None
    missed = 0
    stop_reason = 'max_epochs'
    stopped = False
    training_start = time.monotonic()
    step = 0
    adversary_steps = 0
    autoencoder_steps = 0

    def progress(phase,**extra):
        write_json(directory/'progress.json',{'phase':phase,'training_seconds':time.monotonic()-training_start,**extra})

    def snapshot(epoch):
        return {'identity':request['identity'],'config_sha256':config_digest,'config':config,
                'model_state':model.state_dict(),'optimizer_states':[opt.state_dict() for opt in optimizers],
                'scheduler_states':[scheduler.state_dict() for scheduler in schedulers],
                'epoch':epoch,'best_validation_loss':best,'training_seconds':time.monotonic()-training_start,
                'rng_state':torch.get_rng_state()}

    def stop(signum, frame):
        nonlocal stopped
        stopped = True

    def interrupt(epoch):
        save_checkpoint(directory/'interrupted.pt',snapshot(epoch))
        progress('stopped',stop_reason='signal_or_training_budget')
        raise RuntimeError('chemCPA interrupted; saved interrupted.pt, not a completed prediction')

    signal.signal(signal.SIGTERM,stop)
    signal.signal(signal.SIGINT,stop)
    progress('training',epoch=0)
    for epoch in range(config['max_epochs']):
        model.train()
        sums = np.zeros(3,dtype=float)
        count = 0
        penalties = []
        for indices in minibatches(len(inputs.cases['train']),config['batch_size'],rng,training=True):
            if stopped or time.monotonic()-training_start >= config['max_training_seconds']:
                interrupt(epoch)
            (baseline,fp,doses,expression),drug_ids = context(inputs,'train',indices,device)
            for optimizer in optimizers:
                optimizer.zero_grad(set_to_none=True)
            # Native training reconstructs observed training expression; controls
            # are encoded only for counterfactual prediction, as in author code.
            mean,variance,latent = forward(model,expression,drug_ids,doses)
            reconstruction,adversary,ae_loss,predictions = losses(model,mean,variance,latent,expression,drug_ids)
            if step % hp['adversary_steps'] == 0:
                gradient = torch.autograd.grad(predictions.sum(),latent,create_graph=True)[0]
                penalty = gradient.square().mean()
                objective = adversary+hp['penalty_adversary']*penalty
                if not torch.isfinite(objective):
                    raise ValueError('nonfinite native adversary objective')
                objective.backward()
                optimizers[1].step()
                penalties.append(float(penalty.detach()))
                adversary_steps += 1
            else:
                if not torch.isfinite(ae_loss):
                    raise ValueError('nonfinite native autoencoder objective')
                ae_loss.backward()
                torch.nn.utils.clip_grad_norm_(ae_parameters,1.)
                torch.nn.utils.clip_grad_norm_(model.dosers.parameters(),1.)
                optimizers[0].step()
                optimizers[2].step()
                autoencoder_steps += 1
            sums += np.array([float(reconstruction.detach()),float(adversary.detach()),float(ae_loss.detach())])*len(indices)
            count += len(indices)
            step += 1
            progress('training',epoch=epoch+1,adversary_steps=adversary_steps,autoencoder_steps=autoencoder_steps)
        for scheduler in schedulers:
            scheduler.step()
        model.eval()
        validation_sums = np.zeros(3,dtype=float)
        validation_count = 0
        with torch.no_grad():
            for indices in minibatches(len(inputs.cases['validation']),config['batch_size']):
                if stopped:
                    break
                (_,_,doses,expression),drug_ids = context(inputs,'validation',indices,device)
                mean,variance,latent = forward(model,expression,drug_ids,doses)
                reconstruction,adversary,ae_loss,_ = losses(model,mean,variance,latent,expression,drug_ids)
                validation_sums += np.array([float(reconstruction),float(adversary),float(ae_loss)])*len(indices)
                validation_count += len(indices)
        if stopped:
            interrupt(epoch)
        validation = validation_sums/validation_count
        if not np.isfinite(validation).all():
            raise ValueError('nonfinite native validation objective')
        history.append({'epoch':epoch+1,'reconstruction_loss':float(sums[0]/count),
                        'drug_adversary_BCE':float(sums[1]/count),'autoencoder_objective':float(sums[2]/count),
                        'adversary_gradient_penalty':float(np.mean(penalties)) if penalties else 0.,
                        'validation_reconstruction_loss':float(validation[0]),'validation_drug_adversary_BCE':float(validation[1]),
                        'validation_autoencoder_objective':float(validation[2])})
        if validation[2] < best:
            best = float(validation[2])
            best_epoch = epoch+1
            missed = 0
            save_checkpoint(directory/'best.pt',snapshot(epoch+1))
        else:
            missed += 1
        progress('training',epoch=epoch+1,best_epoch=best_epoch,best_validation_loss=best)
        if missed >= config['early_stop_patience']:
            stop_reason = 'patience'
            break
    training_seconds = time.monotonic()-training_start
    inference_start = time.monotonic()
    write_json(directory/'progress.json',{'phase':'inference','training_seconds':training_seconds})
    selected = torch.load(directory/'best.pt',map_location=device,weights_only=False)
    if selected['identity'] != request['identity'] or selected['config_sha256'] != config_digest:
        raise ValueError('stale chemCPA checkpoint identity')
    model.load_state_dict(selected['model_state'])
    model.eval()
    probe_indices = np.arange(min(2,len(inputs.cases['test'])))
    (baseline,_,doses),drug_ids = context(inputs,'test',probe_indices,device)
    with torch.no_grad():
        mean,_,_ = forward(model,baseline,drug_ids,doses)
        alternate,_,_ = forward(model,torch.zeros_like(baseline),drug_ids,doses)
        dose_alternate,_,_ = forward(model,baseline,drug_ids,doses*2+1.)
        behavior = {'baseline_sensitivity_l2':float(torch.linalg.vector_norm(mean-alternate)),
                    'dose_sensitivity_l2':float(torch.linalg.vector_norm(mean-dose_alternate)),
                    'adapter_parameter_change_l2':float(torch.linalg.vector_norm(vector([model.drug_embeddings.linear])-initial['adapter'])),
                    'predictor_parameter_change_l2':float(torch.linalg.vector_norm(vector([model.encoder,model.decoder,model.drug_embedding_encoder])-initial['predictor'])),
                    'adversary_parameter_change_l2':float(torch.linalg.vector_norm(vector([model.adversary_drugs])-initial['adversary'])),
                    'doser_parameter_change_l2':float(torch.linalg.vector_norm(vector([model.dosers])-initial['doser']))}
        reloaded = make_model(inputs,config,device,seed)
        reloaded.load_state_dict(torch.load(directory/'best.pt',map_location=device,weights_only=False)['model_state'])
        reloaded.eval()
        replayed,_,_ = forward(reloaded,baseline,drug_ids,doses)
        behavior['reload_max_abs_error'] = float(torch.max(torch.abs(mean-replayed)))
    if any(value <= 0 for key,value in behavior.items() if key != 'reload_max_abs_error') or behavior['reload_max_abs_error'] != 0:
        raise ValueError('chemCPA chemical/baseline/dose/adversary/reload behavioral qualification failed')
    del model
    prediction = np.lib.format.open_memmap(directory/'prediction.npy.partial',mode='w+',dtype='float32',
                                          shape=(len(inputs.cases['test']),request['n_genes']))
    with torch.no_grad():
        for indices in minibatches(len(inputs.cases['test']),config['batch_size']):
            if stopped or time.monotonic()-inference_start >= config['inference_max_seconds']:
                raise RuntimeError('chemCPA inference interrupted/budget exceeded; no complete prediction')
            (baseline,_,doses),drug_ids = context(inputs,'test',indices,device)
            mean,_,_ = forward(reloaded,baseline,drug_ids,doses)
            delta = (mean-baseline).cpu().numpy()
            if not np.isfinite(delta).all():
                raise ValueError('nonfinite whole-axis chemCPA prediction')
            prediction[indices] = delta
    prediction.flush()
    del prediction
    (directory/'prediction.npy.partial').replace(directory/'prediction.npy')
    write_json(directory/'case_ids.json',[r['case_id'] for r in inputs.cases['test']])
    inference_seconds = time.monotonic()-inference_start
    result = {'identity':request['identity'],'author_objective':'GaussianNLL - reg_adversary * drug_BCE',
              'best_epoch':best_epoch,'best_validation_loss':best,'history':history,'stop_reason':stop_reason,
              'training_seconds':training_seconds,'inference_seconds':inference_seconds,
              'trainable_parameters':sum(p.numel() for p in reloaded.parameters()),
              'adapter_parameters':sum(p.numel() for p in reloaded.drug_embeddings.linear.parameters()),
              'adversary_steps':adversary_steps,'autoencoder_steps':autoencoder_steps,'behavior':behavior,
              'blobs':{'best.pt':sha256(directory/'best.pt')},'config':config,'config_sha256':config_digest,
              'selection':'minimum native combined autoencoder objective on validation, never test response metrics',
              'training_method':'native treated-expression reconstruction plus adversarial disentanglement',
              'covariates':'none: frozen tenv request has no variable cell-type covariate',
              'output':'native expression mean minus approved same-plate control baseline',
              'environment':{'python':sys.version,'executable':sys.executable,'numpy':np.__version__,'torch':torch.__version__,
                             'cuda':torch.version.cuda,'device':str(device)},
              'peak_rss_kbytes':int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
              'peak_gpu_bytes':int(torch.cuda.max_memory_allocated()) if torch.cuda.is_available() else 0,
              'peak_gpu_reserved_bytes':int(torch.cuda.max_memory_reserved()) if torch.cuda.is_available() else 0}
    write_json(directory/'checkpoint.json',result)
    write_json(directory/'progress.json',{'phase':'complete','training_seconds':training_seconds,'inference_seconds':inference_seconds})
    print(json.dumps({key:result[key] for key in ('best_epoch','best_validation_loss','training_seconds','inference_seconds','behavior')}),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request',required=True)
    parser.add_argument('--config',default=str(Path(__file__).resolve().parents[1]/'models/chemcpa/config.json'))
    args = parser.parse_args()
    fit_predict(read_json(args.request),read_json(args.config))
