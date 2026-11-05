"""Separate native PerturbDiff synthetic setup proof; no Morgan adapter or tenv data."""
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import tempfile
import time
sys.dont_write_bytecode=True
for name in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(name,'2')
import numpy as np
import torch
from omegaconf import OmegaConf
REPO=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
UPSTREAM=REPO/'pert_model/PerturbDiff/source'
sys.path.insert(0,str(UPSTREAM))
from src.models.lightning.lightning_module import PlModel
from src.models.lightning.lightning_factories import get_optimizer


def main():
    torch.set_num_threads(2); torch.manual_seed(17); np.random.seed(17)
    device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    if device.type=='cuda':
        torch.cuda.reset_peak_memory_stats()
    start=time.monotonic()
    model_cfg=OmegaConf.load(UPSTREAM/'configs/model/base.yaml')
    model_cfg.class_emb_gather_strategy='concat'; model_cfg.class_emb_hidden_dimension=128
    model_cfg.dataset_dict={'tahoe_native_synthetic':0}; model_cfg.p_drop_control=0.
    cov_cfg=OmegaConf.create({'pert_encoding':'onehot','drug_encoding':'onehot','replogle_gene_encoding':'onehot',
                             'hidden_dim':128,'output_dim':2000,'num_pert':2,'celltype_encoding':'onehot',
                             'num_celltype':1,'batch_encoding':None})
    optim=OmegaConf.create({'micro_batch_size':2,'blur':.05,'use_mse_loss':True,'MMD_loss_factor':1.,
                            'optimizer':{'lr':.0002,'adam_w_mode':True,'eps':1e-8,'betas':[.9,.98],'weight_decay':.01}})
    model=PlModel(cov_cfg,model_cfg,None,optim,OmegaConf.create({'log_every_n_steps':20}),['validation']).to(device)
    optimizer=get_optimizer([{'params':model.parameters(),'weight_decay':.01}],optim.optimizer)
    def batch():
        control=torch.rand(2,3,2000,device=device)*.1
        return {'cont_emb':control,'pert_emb':control+torch.rand_like(control)*.02,
                'cov_pert':torch.tensor([[0],[1]],device=device),'cov_celltype':torch.zeros(2,1,dtype=torch.long,device=device),
                'cov_batch':torch.zeros(2,1,dtype=torch.long,device=device),'ds_name':['tahoe_native_synthetic']*2}
    losses=[]
    for step in range(2):
        optimizer.zero_grad(set_to_none=True)
        result=model._compute_loss(batch())
        loss=(result['loss']*result['weights']).mean()
        if not torch.isfinite(loss):
            raise ValueError('native setup nonfinite objective')
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),1.); optimizer.step()
        losses.append(float(loss.detach()))
    model.eval(); probe=batch()
    cond={'batch_emb':model._encode_covariates(probe),'cont_emb':probe['cont_emb'],'gene_emb':None,
          'ds_name':probe['ds_name']}
    noise=torch.zeros_like(probe['cont_emb'])
    def sample():
        return model.diffusion.ddim_sample_loop(model.model,tuple(noise.shape),self_condition=cond,
                 noise=noise,device=device,start_time=100,eta=0.,guidance_strength=1.,clip_denoised=True)[0]
    with torch.no_grad(),tempfile.TemporaryDirectory() as temp:
        before=sample(); path=Path(temp)/'native.pt'; torch.save(model.state_dict(),path)
        checkpoint_sha=hashlib.sha256(path.read_bytes()).hexdigest()
        model.load_state_dict(torch.load(path,map_location=device,weights_only=True))
        after=sample()
    error=float((before-after).abs().max())
    if error!=0 or not torch.isfinite(after).all():
        raise ValueError('native reload/sampling failed')
    result={'purpose':'native setup only; synthetic unmodified model, no Morgan or tenv qualification',
            'source_commit':'31be89cb8cce5dddef0782bac29d5c748ce55931','shape':list(after.shape),
            'losses':losses,'objective':'native GaussianDiffusion MSE + energy MMD','diffusion_steps':1000,
            'ddim_start_time':100,'reload_max_abs_error':error,'checkpoint_sha256':checkpoint_sha,
            'trainable_parameters':sum(p.numel() for p in model.parameters()),'seconds':time.monotonic()-start,
            'peak_rss_kbytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'peak_gpu_bytes':torch.cuda.max_memory_allocated() if device.type=='cuda' else 0,
            'environment':{'executable':sys.executable,'python':sys.version,'torch':torch.__version__,
                           'numpy':np.__version__,'device':str(device)}}
    output=Path(__file__).resolve().parent/'output'; output.mkdir(exist_ok=True)
    (output/'meta.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
