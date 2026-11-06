"""Isolated CMonge v2: fixed full validation and predeclared stopping."""
import argparse
import os
from pathlib import Path
import signal
import sys
sys.dont_write_bytecode = True
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
from cmonge_v2_protocol import PROTOCOL, EVAL_EVERY, Selection, remaining_seconds

class AdaptedTransport(nn.Module):
    genes: int
    @nn.compact
    def __call__(self, baseline, morgan, logdose):
        chemistry=nn.Dense(50,name='morgan_adapter')(morgan)
        context=jnp.concatenate((chemistry,logdose),axis=1)
        return ConditionalPerturbationNetwork(dim_hidden=(64,64,64,64),dim_data=self.genes,
                   dim_cond=51,dim_cond_map=50,context_entity_bonds=((0,50),(0,51)),name='native_transport')(baseline,context,2)


def run(request, config):
    remaining=remaining_seconds(request)
    output=Path(request['output_directory']); files=request['predictor_files']
    if config.get('protocol')!=PROTOCOL:raise ValueError('explicit CMonge v2 protocol required')
    qualification=bool(config.get('qualification',False))
    n_genes=int(request['n_genes'])
    if not qualification and n_genes!=18164: raise ValueError('formal CMonge requires full gene axis')
    seed=request['identity']['seed']; rng=np.random.default_rng(seed)
    validation_rng=np.random.default_rng(seed+7927)
    iterations=int(config.get('iterations',10000)); eval_every=EVAL_EVERY
    if int(config.get('eval_every',EVAL_EVERY))!=EVAL_EVERY:raise ValueError('v2 cadence fixed at50')
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
    if not qualification and len(validation[2])!=138:raise ValueError('full138 validation groups required')
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
    signal.signal(signal.SIGINT,lambda *_:stopped.__setitem__(0,True))
    started=time.monotonic(); deadline=started+remaining; best=None; best_loss=float('inf'); history=[]
    selection=Selection(); stop_reason='maximum_iterations'
    def guard():
        if stopped[0] or time.monotonic()>=deadline:
            write_json(output/'progress.json',{'phase':'training','training_seconds':time.monotonic()-started,
                       'stop_reason':'signal' if stopped[0] else 'cumulative_budget_exhausted'})
            raise SystemExit(143 if stopped[0] else 124)
    write_json(output/'progress.json',{'phase':'training','training_seconds':0.})
    for iteration in range(iterations):
        guard()
        batch_data=sample(training,int(rng.integers(len(training[2]))),rng)
        state,total,parts=step(state,batch_data)
        total=float(total)
        if not np.isfinite(total): raise ValueError('nonfinite native author training objective')
        write_json(output/'progress.json',{'phase':'training','training_seconds':time.monotonic()-started,'iteration':iteration+1,'author_training_loss':total})
        if (iteration+1)%eval_every==0 or iteration+1==iterations:
            # Reset only the validation stream: identical draws for every group
            # at every event. Training RNG and native batch objective are unchanged.
            fixed_validation_rng=np.random.default_rng(seed+7927)
            values=[]; validation_started=time.monotonic()
            for group_index in range(len(validation[2])):
                guard()
                group_batch=sample(validation,group_index,fixed_validation_rng)
                value=float(loss_eval(state.params,*group_batch)[0])
                if not np.isfinite(value):raise ValueError('nonfinite validation objective')
                values.append(value)
                del group_batch
            guard()
            value=float(np.mean(values,dtype=np.float64))
            improved,early_stop=selection.observe(iteration+1,value,final=iteration+1==iterations)
            history.append({'iteration':iteration+1,'validation_loss':value,
                'validation_conditions_evaluated':len(values),'training_loss':total,
                'fitting_loss':float(parts[0]),'regularizer':float(parts[1]),
                'validation_seconds':time.monotonic()-validation_started,
                'early_stop_bad_checks':selection.bad_checks})
            # Durable event artifacts survive timeout; these snapshots do not
            # claim exact continuation (RNG/scheduler state is not restorable).
            (output/'latest.msgpack').write_bytes(serialization.to_bytes(state.params))
            if improved:
                best_loss=value;best=jax.tree_util.tree_map(lambda x:np.array(x),state.params)
                (output/'best.msgpack').write_bytes(serialization.to_bytes(best))
            write_json(output/'history.json',{'protocol':PROTOCOL,'identity':request['identity'],'history':history,
                       'snapshot_semantics':'prediction-only-not-resumable'})
            guard()
            if early_stop:
                stop_reason='validation_patience';break
    guard()
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
    if not np.isfinite([adapter_change,baseline_sensitivity,dose_sensitivity,reload_error]).all() or not adapter_change>0 or reload_error!=0:
        raise ValueError('nonfinite diagnostics or adapter/reload behavior check failed')
    # Admission proves that inputs participate. A fitted formal checkpoint may
    # be locally insensitive on this one condition; retain that model behavior
    # as a diagnostic instead of censoring an otherwise valid prediction run.
    if qualification and (not baseline_sensitivity>0 or not dose_sensitivity>0):
        raise ValueError('qualification baseline/dose behavior check failed')
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
    write_json(output/'checkpoint.json',{'identity':request['identity'],'blobs':{'best.msgpack':sha256(output/'best.msgpack'),'latest.msgpack':sha256(output/'latest.msgpack'),'history.json':sha256(output/'history.json')},
               'training_seconds':training_seconds,'inference_seconds':time.monotonic()-inference_start,
               'author_objective':'Sinkhorn divergence + Monge gap (Euclidean)','best_validation_loss':best_loss,
               'selection':'minimum equal-group full138 fixed-draw native author objective','history':history,
               'protocol':PROTOCOL,'stop_reason':stop_reason,'best_iteration':selection.best_iteration,
               'snapshot_semantics':'prediction-only-not-resumable','prior_training_seconds':request['prior_training_seconds'],
               'early_stopping':{'minimum_iterations':500,'patience':10,'relative_min_delta':.001,'earliest_stop_iteration':1000},
               'adapter_parameter_change':adapter_change,'baseline_sensitivity':baseline_sensitivity,'dose_sensitivity':dose_sensitivity,
               'input_sensitivity_check':'qualification admission gate' if qualification else 'formal single-condition diagnostic; zero retained',
               'reload_max_absolute_error':reload_error,'trainable_parameters':sum(int(x.size) for x in jax.tree_util.tree_leaves(restored)),
               'adapter_parameters':2048*50+50,'adapter_dimensions':[2048,50],'batch_size':batch,'iterations':iteration+1,'maximum_iterations':iterations,
               'validation_groups':len(validation[2]),'validation_batch_policy':'all groups sequentially, equal group mean, identical seeded within-group draws each event',
               'qualification_only':qualification,'gene_count':n_genes,'backend':jax.default_backend(),
               'device_kind':jax.devices()[0].device_kind,'device_memory_stats':jax.devices()[0].memory_stats(),
               'recipe_deviations':['Shared Morgan linear map replaces precomputed RDKit descriptor representation.',
                    'Direct native transport on complete gene space; optional separately trained AE omitted.',
                    'Same-plate DMSO mean replaces distribution of individual control cells.',
                    'Explicit log(float micromolar dose) preserves fractional concentrations.',
                    'V2 selects minimum fixed-draw full-validation native objective every50 iterations; patience10 after500 iterations with relative0.001 improvement; final iteration also validated.']})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--request',required=True);parser.add_argument('--config',required=True)
    args=parser.parse_args();run(read_json(args.request),read_json(args.config))
