"""CMonge native conditional transport with jointly trained Morgan linear map."""
import argparse
import os
from pathlib import Path
import signal
import sys
import time
os.environ.setdefault('JAX_PLATFORM_NAME','cpu')
os.environ.setdefault('XLA_PYTHON_CLIENT_PREALLOCATE','false')
os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2')
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
sys.path.insert(0,str(ROOT/'pert_model/CMonge/source'))
import jax
import jax.numpy as jnp
import flax.linen as nn
from flax import serialization
from flax.training.train_state import TrainState
import numpy as np
import optax
from cmonge.models.nn import ConditionalPerturbationNetwork
from cmonge.metrics import fitting_loss, regularizer
from artifacts import read_json, write_json, rows, sha256

class AdaptedTransport(nn.Module):
    genes: int
    @nn.compact
    def __call__(self, baseline, morgan, logdose):
        chemistry=nn.Dense(50,name='morgan_adapter')(morgan)
        context=jnp.concatenate((chemistry,logdose),axis=1)
        return ConditionalPerturbationNetwork(dim_hidden=(64,64,64,64),dim_data=self.genes,
                   dim_cond=51,dim_cond_map=50,context_entity_bonds=((0,50),(0,51)),name='native_transport')(baseline,context,2)


def run(request, config):
    output=Path(request['output_directory']); files=request['predictor_files']
    qualification=bool(config.get('qualification',False))
    n_genes=int(request['n_genes'])
    if not qualification and n_genes!=18164: raise ValueError('formal CMonge requires full gene axis')
    seed=request['identity']['seed']; rng=np.random.default_rng(seed)
    iterations=int(config.get('iterations',10000)); eval_every=int(config.get('eval_every',100))
    batch=max(1,int(config.get('batch_size',512)*float(os.environ.get('TENV_BATCH_SCALE','1'))))
    if iterations<1 or batch<1 or eval_every<1: raise ValueError('invalid training budget')
    controls=np.load(files['controls'],mmap_mode='r'); morgan=np.load(files['morgan'],mmap_mode='r')
    drug_index={row['drug_key']:int(row['morgan_index']) for row in rows(files['drugs'])}
    def data(role):
        cases=rows(files[f'{role}_cases']); expression=np.load(files[f'{role}_expression'],mmap_mode='r')
        groups=sorted({r['group_id'] for r in cases})
        if qualification: groups=groups[:int(config.get('max_groups',3))]
        return cases,expression,[[i for i,r in enumerate(cases) if r['group_id']==group] for group in groups]
    training=data('train'); validation=data('validation')
    def sample(dataset, group, generator):
        cases,delta,groups=dataset
        indices=generator.choice(groups[group],size=batch,replace=True)
        selected=[cases[i] for i in indices]
        dose=np.array([float(r['dose_um']) for r in selected],dtype='float32')
        if np.any(dose<=0): raise ValueError('CMonge log-dose requires positive actual doses')
        baseline=np.array([controls[int(r['control_index'])] for r in selected],dtype='float32')
        chemical=np.array([morgan[drug_index[r['drug_key']]] for r in selected],dtype='float32')
        target=baseline+np.asarray(delta[indices],dtype='float32')
        return tuple(jnp.asarray(v) for v in (baseline,chemical,np.log(dose)[:,None],target))
    model=AdaptedTransport(n_genes)
    first=sample(training,0,rng)
    params=model.init(jax.random.PRNGKey(seed),*first[:3])['params']
    original_adapter=np.asarray(params['morgan_adapter']['kernel']).copy()
    schedule=optax.cosine_decay_schedule(1e-4,iterations,alpha=1e-2)
    state=TrainState.create(apply_fn=model.apply,params=params,tx=optax.adamw(schedule,weight_decay=1e-5))
    def loss(params, baseline, chemical, dose, target):
        mapped=model.apply({'params':params},baseline,chemical,dose)
        fit=fitting_loss(target,mapped,epsilon_fitting=1.)
        reg=regularizer(baseline,mapped,epsilon_regularizer=.01,cost='euclidean')
        return fit+reg,(fit,reg)
    @jax.jit
    def step(state,batch):
        (total,parts),grads=jax.value_and_grad(loss,has_aux=True)(state.params,*batch)
        return state.apply_gradients(grads=grads),total,parts
    loss_eval=jax.jit(loss)
    predict=jax.jit(lambda p,x,c,d:model.apply({'params':p},x,c,d))
    stopped=[False]
    signal.signal(signal.SIGTERM,lambda *_:stopped.__setitem__(0,True))
    started=time.monotonic(); best=None; best_loss=float('inf'); history=[]
    write_json(output/'progress.json',{'phase':'training','training_seconds':0.})
    for iteration in range(iterations):
        if stopped[0]:
            (output/'interrupted.msgpack').write_bytes(serialization.to_bytes(state))
            write_json(output/'progress.json',{'phase':'training','training_seconds':time.monotonic()-started,'stop_reason':'SIGTERM checkpoint saved'})
            raise SystemExit(143)
        batch_data=sample(training,int(rng.integers(len(training[2]))),rng)
        state,total,parts=step(state,batch_data)
        total=float(total)
        if not np.isfinite(total): raise ValueError('nonfinite native author training objective')
        write_json(output/'progress.json',{'phase':'training','training_seconds':time.monotonic()-started,'iteration':iteration+1,'author_training_loss':total})
        if (iteration+1)%eval_every==0 or iteration+1==iterations:
            # Regenerate one fixed group batch at a time. Caching all validation
            # batches at 18k genes would exceed the formal process memory cap.
            losses=[]
            for group_index in range(len(validation[2])):
                fixed_rng=np.random.default_rng(seed+7927+group_index*104729)
                group_batch=sample(validation,group_index,fixed_rng)
                losses.append(float(loss_eval(state.params,*group_batch)[0]))
                del group_batch
            value=float(np.mean(losses))
            if not np.isfinite(value): raise ValueError('nonfinite native author validation objective')
            history.append({'iteration':iteration+1,'validation_loss':value,'training_loss':total,'fitting_loss':float(parts[0]),'regularizer':float(parts[1])})
            if value<best_loss:
                best_loss=value;best=jax.tree_util.tree_map(lambda x:np.array(x),state.params)
                (output/'best.msgpack').write_bytes(serialization.to_bytes(best))
    training_seconds=time.monotonic()-started
    write_json(output/'progress.json',{'phase':'inference','training_seconds':training_seconds})
    inference_start=time.monotonic()
    if best is None:raise ValueError('no author validation checkpoint')
    before=np.asarray(predict(best,*first[:3]))
    restored=serialization.from_bytes(params,(output/'best.msgpack').read_bytes())
    after=np.asarray(predict(restored,*first[:3]))
    reload_error=float(np.max(np.abs(before-after)))
    adapter_change=float(np.linalg.norm(np.asarray(restored['morgan_adapter']['kernel'])-original_adapter))
    baseline_sensitivity=float(np.max(np.abs(np.asarray(predict(restored,first[0]+.1,*first[1:3]))-after)))
    dose_sensitivity=float(np.max(np.abs(np.asarray(predict(restored,*first[:2],first[2]+1.))-after)))
    if not adapter_change>0 or not baseline_sensitivity>0 or not dose_sensitivity>0 or reload_error!=0:
        raise ValueError('adapter/baseline/dose/reload behavior check failed')
    cases=rows(files['test_cases'])
    prediction=np.lib.format.open_memmap(output/'prediction.npy',mode='w+',dtype='float32',shape=(len(cases),n_genes))
    for start in range(0,len(cases),32):
        selected=cases[start:start+32]
        baseline=np.array([controls[int(r['control_index'])] for r in selected],dtype='float32')
        chemical=np.array([morgan[drug_index[r['drug_key']]] for r in selected],dtype='float32')
        dose=np.array([float(r['dose_um']) for r in selected],dtype='float32')
        if np.any(dose<=0):raise ValueError('nonpositive test dose')
        prediction[start:start+len(selected)]=np.asarray(predict(restored,jnp.asarray(baseline),jnp.asarray(chemical),jnp.asarray(np.log(dose)[:,None])))-baseline
    prediction.flush()
    write_json(output/'case_ids.json',[r['case_id'] for r in cases])
    write_json(output/'checkpoint.json',{'identity':request['identity'],'blobs':{'best.msgpack':sha256(output/'best.msgpack')},
               'training_seconds':training_seconds,'inference_seconds':time.monotonic()-inference_start,
               'author_objective':'Sinkhorn divergence + Monge gap (Euclidean)','best_validation_loss':best_loss,
               'selection':'minimum mean held-out-group native author objective','history':history,
               'adapter_parameter_change':adapter_change,'baseline_sensitivity':baseline_sensitivity,'dose_sensitivity':dose_sensitivity,
               'reload_max_absolute_error':reload_error,'trainable_parameters':sum(int(x.size) for x in jax.tree_util.tree_leaves(restored)),
               'adapter_parameters':2048*50+50,'adapter_dimensions':[2048,50],'batch_size':batch,'iterations':iterations,
               'validation_groups':len(validation[2]),'validation_batch_policy':'one streamed group; fixed per-group seed',
               'qualification_only':qualification,'gene_count':n_genes,'backend':jax.default_backend(),
               'recipe_deviations':['Shared Morgan linear map replaces precomputed RDKit descriptor representation.',
                    'Direct native transport on complete gene space; optional separately trained AE omitted.',
                    'Same-plate DMSO mean replaces distribution of individual control cells.',
                    'Explicit log(float micromolar dose) preserves fractional concentrations.',
                    'Checkpoint selection evaluates native total validation objective across fixed held-out groups.']})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--request',required=True);parser.add_argument('--config',required=True)
    args=parser.parse_args();run(read_json(args.request),read_json(args.config))
