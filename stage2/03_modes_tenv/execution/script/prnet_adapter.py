"""Isolated native PRnet PGM with the agreed joint Morgan linear adapter."""
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
    os.environ.setdefault(_name, '2')

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F
from artifacts import read_json, write_json, rows, sha256

def repo_root():
    for candidate in Path(__file__).resolve().parents:
        if (candidate/'.git').exists():
            return candidate
    raise RuntimeError('repository root not found')

REPO = repo_root()
UPSTREAM = REPO/'pert_model/PRnet/source'
sys.path.insert(0, str(UPSTREAM))
from models.PRnet import PGM


class JointPRnet(nn.Module):
    def __init__(self, n_genes, config):
        super().__init__()
        self.adapter = nn.Linear(2048, config['drug_dimension'])
        self.predictor = PGM(n_genes, config['comb_dimension'], config['noise_dimension'],
                             config['hidden_layer_sizes'], config['z_dimension'],
                             config['adaptor_layer_sizes'], config['drug_dimension'],
                             config['dropout'])
        # Upstream PRnetTrainer.weight_init overrides native BatchNorm defaults.
        for module in self.predictor.modules():
            if isinstance(module,nn.BatchNorm1d):
                nn.init.normal_(module.weight,1.,.02)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)

    def forward(self, baseline, fingerprints, dose, noise):
        # Author Drug_dose_encoder scales chemistry before native CombAdaptor.
        chemical = self.adapter(fingerprints) * torch.log10(dose[:, None] + 1.)
        return self.predictor(baseline, chemical, noise)


from predictor_inputs import PredictorInputs, minibatches


def parameter_vector(module):
    return torch.cat([p.detach().cpu().reshape(-1) for p in module.parameters()])


def save_torch(path, value):
    temporary = path.with_suffix(path.suffix+'.partial')
    torch.save(value,temporary)
    temporary.replace(path)


def fit_predict(request, config):
    config = dict(config)
    config['batch_size'] = max(2, int(config['batch_size'] * float(os.environ.get('TENV_BATCH_SCALE', '1'))))
    if config['max_epochs'] < 1 or config['batch_size'] < 2 or config['threads'] not in (1,2):
        raise ValueError('invalid bounded PRnet recipe')
    torch.set_num_threads(config['threads'])
    seed = int(request['identity']['seed'])
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.cuda.reset_peak_memory_stats()
    device = torch.device(('cuda' if torch.cuda.is_available() else 'cpu') if config['device']=='auto' else config['device'])
    directory = Path(request['output_directory'])
    directory.mkdir(parents=True,exist_ok=True)
    inputs = PredictorInputs(request)
    model = JointPRnet(request['n_genes'],config).to(device)
    initial_adapter = parameter_vector(model.adapter)
    initial_predictor = parameter_vector(model.predictor)
    optimizer = torch.optim.Adam(model.parameters(),lr=config['learning_rate'],weight_decay=config['weight_decay'])
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer,mode='min',factor=config['scheduler_factor'],
                    patience=config['scheduler_patience'],min_lr=config['minimum_learning_rate'])
    objective = nn.GaussianNLLLoss()
    rng = np.random.default_rng(seed)
    history = []
    training_start = time.monotonic()
    best = float('inf')
    best_epoch = None
    missed = 0
    stopped = False
    stop_reason = 'max_epochs'
    config_digest = hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()

    def progress(phase, **extra):
        elapsed = time.monotonic()-training_start
        write_json(directory/'progress.json',{'phase':phase,'training_seconds':elapsed,**extra})

    def snapshot(epoch):
        return {'identity':request['identity'],'config':config,'config_sha256':config_digest,
                'model_state':model.state_dict(),'optimizer_state':optimizer.state_dict(),
                'scheduler_state':scheduler.state_dict(),'epoch':epoch,'best_validation_loss':best,
                'training_seconds':time.monotonic()-training_start,'rng_state':torch.get_rng_state()}

    def request_stop(signum, frame):
        nonlocal stopped
        stopped = True

    signal.signal(signal.SIGTERM,request_stop)
    signal.signal(signal.SIGINT,request_stop)
    progress('training',epoch=0)
    for epoch in range(config['max_epochs']):
        model.train()
        total = 0.
        count = 0
        for indices in minibatches(len(inputs.cases['train']),config['batch_size'],rng,training=True):
            if stopped or time.monotonic()-training_start >= config['max_training_seconds']:
                stopped = True
                break
            baseline,fingerprint,dose,target = inputs.batch('train',indices,device)
            optimizer.zero_grad(set_to_none=True)
            output = model(baseline,fingerprint,dose,torch.randn(len(indices),config['noise_dimension'],device=device))
            mean,variance = output.chunk(2,dim=1)
            loss = objective(mean,target,F.softplus(variance))
            if not torch.isfinite(loss):
                raise ValueError('nonfinite author training objective')
            loss.backward()
            optimizer.step()
            total += float(loss.detach())*len(indices)
            count += len(indices)
            progress('training',epoch=epoch+1)
        if stopped:
            save_torch(directory/'interrupted.pt',snapshot(epoch))
            progress('stopped',stop_reason='signal_or_training_budget')
            raise RuntimeError('PRnet interrupted: saved interrupted.pt; no completed prediction')
        model.eval()
        validation_total = 0.
        validation_count = 0
        # Repeat a fixed validation noise stream each epoch to compare checkpoints.
        valid_generator = torch.Generator(device=device).manual_seed(seed+100000)
        with torch.no_grad():
            for indices in minibatches(len(inputs.cases['validation']),config['batch_size']):
                if stopped:
                    break
                baseline,fingerprint,dose,target = inputs.batch('validation',indices,device)
                noise = torch.randn(len(indices),config['noise_dimension'],generator=valid_generator,device=device)
                mean,variance = model(baseline,fingerprint,dose,noise).chunk(2,dim=1)
                loss = objective(mean,target,F.softplus(variance))
                validation_total += float(loss)*len(indices)
                validation_count += len(indices)
        if stopped:
            save_torch(directory/'interrupted.pt',snapshot(epoch))
            progress('stopped',stop_reason='signal')
            raise RuntimeError('PRnet interrupted during validation: saved interrupted.pt')
        validation = validation_total/validation_count
        if not np.isfinite(validation):
            raise ValueError('nonfinite author validation objective')
        scheduler.step(validation)
        history.append({'epoch':epoch+1,'training_loss':total/count,'validation_loss':validation,
                        'learning_rate':optimizer.param_groups[0]['lr']})
        if validation < best:
            best = validation
            best_epoch = epoch+1
            missed = 0
            save_torch(directory/'best.pt',snapshot(epoch+1))
        else:
            missed += 1
        progress('training',epoch=epoch+1,best_epoch=best_epoch,best_validation_loss=best)
        if missed >= config['early_stop_patience']:
            stop_reason = 'author_patience'
            break
    training_seconds = time.monotonic()-training_start
    inference_start = time.monotonic()
    write_json(directory/'progress.json',{'phase':'inference','training_seconds':training_seconds})
    # Validate the persisted selected state and compare a real forward after reload.
    checkpoint = torch.load(directory/'best.pt',map_location=device,weights_only=False)
    if checkpoint['identity'] != request['identity'] or checkpoint['config_sha256'] != config_digest:
        raise ValueError('stale selected checkpoint identity')
    model.load_state_dict(checkpoint['model_state'])
    model.eval()
    probe = inputs.batch('test',np.arange(min(2,len(inputs.cases['test']))),device)
    probe_noise = torch.zeros(len(probe[0]),config['noise_dimension'],device=device)
    with torch.no_grad():
        original = model(*probe,probe_noise).detach()
        changed_baseline = model(torch.zeros_like(probe[0]),probe[1],probe[2],probe_noise)
        baseline_sensitivity = float(torch.linalg.vector_norm(
            original[:,:request['n_genes']]-changed_baseline[:,:request['n_genes']]))
        adapter_change = float(torch.linalg.vector_norm(parameter_vector(model.adapter)-initial_adapter))
        predictor_change = float(torch.linalg.vector_norm(parameter_vector(model.predictor)-initial_predictor))
        reloaded = JointPRnet(request['n_genes'],config).to(device)
        reloaded.load_state_dict(torch.load(directory/'best.pt',map_location=device,weights_only=False)['model_state'])
        reloaded.eval()
        replayed = reloaded(*probe,probe_noise)
        reload_error = float(torch.max(torch.abs(original-replayed)))
    if not adapter_change > 0 or not predictor_change > 0 or not baseline_sensitivity > 0 or reload_error != 0:
        raise ValueError('PRnet adapter/baseline/reload behavioral gate failed')
    del model
    prediction = np.lib.format.open_memmap(directory/'prediction.npy.partial',mode='w+',dtype='float32',
                   shape=(len(inputs.cases['test']),request['n_genes']))
    predict_generator = torch.Generator(device=device).manual_seed(seed+200000)
    with torch.no_grad():
        for indices in minibatches(len(inputs.cases['test']),config['batch_size']):
            if stopped or time.monotonic()-inference_start >= config['inference_max_seconds']:
                raise RuntimeError('PRnet inference interrupted/budget exceeded; partial prediction is not completion')
            baseline,fingerprint,dose = inputs.batch('test',indices,device)
            noise = torch.randn(len(indices),config['noise_dimension'],generator=predict_generator,device=device)
            mean = reloaded(baseline,fingerprint,dose,noise)[:,:request['n_genes']]
            delta = (mean-baseline).cpu().numpy()
            if not np.isfinite(delta).all():
                raise ValueError('nonfinite full-axis prediction')
            prediction[indices] = delta
    prediction.flush()
    del prediction
    (directory/'prediction.npy.partial').replace(directory/'prediction.npy')
    write_json(directory/'case_ids.json',[case['case_id'] for case in inputs.cases['test']])
    inference_seconds = time.monotonic()-inference_start
    result = {'identity':request['identity'],'author_objective':'GaussianNLLLoss',
              'best_epoch':best_epoch,'best_validation_loss':best,'history':history,'stop_reason':stop_reason,
              'training_seconds':training_seconds,'inference_seconds':inference_seconds,
              'trainable_parameters':sum(p.numel() for p in reloaded.parameters()),
              'adapter_parameters':sum(p.numel() for p in reloaded.adapter.parameters()),
              'predictor_parameters':sum(p.numel() for p in reloaded.predictor.parameters()),
              'behavior':{'adapter_parameter_change_l2':adapter_change,'predictor_parameter_change_l2':predictor_change,
                          'baseline_sensitivity_l2':baseline_sensitivity,'reload_max_abs_error':reload_error},
              'selection':'minimum author Gaussian NLL on validation; fixed validation noise stream',
              'output':'native conditional Gaussian mean minus approved same-plate baseline',
              'blobs':{'best.pt':sha256(directory/'best.pt')},'config':config,'config_sha256':config_digest,
              'environment':{'python':sys.version,'executable':sys.executable,'torch':torch.__version__,
                             'cuda':torch.version.cuda,'numpy':np.__version__,'device':str(device)},
              'peak_rss_kbytes':int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
              'peak_gpu_bytes':int(torch.cuda.max_memory_allocated()) if torch.cuda.is_available() else 0,
              'peak_gpu_reserved_bytes':int(torch.cuda.max_memory_reserved()) if torch.cuda.is_available() else 0}
    write_json(directory/'checkpoint.json',result)
    write_json(directory/'progress.json',{'phase':'complete','training_seconds':training_seconds,'inference_seconds':inference_seconds})
    print(json.dumps({k:result[k] for k in ('best_epoch','best_validation_loss','stop_reason','training_seconds','inference_seconds','behavior')}),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request',required=True)
    parser.add_argument('--config',default=str(Path(__file__).resolve().parents[1]/'models/prnet/config.json'))
    args = parser.parse_args()
    fit_predict(read_json(args.request),read_json(args.config))
