"""Artifact CLI: isolated predictors, frozen completion barrier, offline evaluation."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from artifacts import (SEEDS, TERMINAL, sha256, read_json, write_json, rows, write_csv,
                       run_directory, run_identity, validate_prediction, validate_run, command_evidence, checkpoint_artifacts, verify_input_files)
from metrics import METRICS, RECALL_GRID, CURVE_CONVENTION, score_group


def predictor_files(input_root, fold):
    root = Path(input_root).resolve(); directory = root / f'fold_{fold}'
    # Explicit allowlist, never root manifest, global well indices or evaluator truth.
    mapping = {'train_cases':directory/'train.csv', 'validation_cases':directory/'validation.csv',
               'test_cases':directory/'test_cases.csv', 'train_expression':directory/'train_expression.npy',
               'validation_expression':directory/'validation_expression.npy',
               'controls':root/'controls.npy', 'morgan':root/'morgan.npy',
               'drugs':root/'drugs.csv', 'genes':root/'genes.csv'}
    allowed = {'case_id','drug_key','dose_um','control_index','group_id'}
    for role in ('train_cases','validation_cases','test_cases'):
        case_rows = rows(mapping[role])
        if not case_rows or set(case_rows[0]) != allowed:
            raise ValueError('predictor case leakage or invalid case schema')
        if len({r['case_id'] for r in case_rows}) != len(case_rows):
            raise ValueError('duplicate predictor case identity')
    drugs = {role:{r['drug_key'] for r in rows(mapping[role])} for role in ('train_cases','validation_cases','test_cases')}
    roles=list(drugs)
    drug_rows=rows(mapping['drugs'])
    block_lookup={r['drug_key']:r.get('parent_connectivity',r['drug_key']) for r in drug_rows}
    for first in range(len(roles)):
        for second in range(first+1,len(roles)):
            left,right=drugs[roles[first]],drugs[roles[second]]
            if left & right:
                raise ValueError('drug leakage across predictor roles')
            if any(drug not in block_lookup for drug in left|right):
                raise ValueError('unknown drug identity')
            if {block_lookup[d] for d in left} & {block_lookup[d] for d in right}:
                raise ValueError('structure-block leakage across predictor roles')
    return {key:str(value) for key,value in mapping.items()}


def run_one(input_root, output_root, model, fold, seed, command=None):
    verify_input_files(input_root)
    directory=run_directory(output_root,model,fold,seed)
    if (directory/'run.json').exists() and read_json(directory/'run.json').get('status')=='completed':
        return validate_run(input_root,output_root,model,fold,seed,command)
    files=predictor_files(input_root,fold)
    identity=run_identity(input_root,model,fold,seed,command)
    request={'identity':identity,'predictor_files':files,'target_kind':'delta',
             'n_genes':len(rows(files['genes'])),'output_directory':str(directory.resolve()),
             'command_evidence':command_evidence(command),
             'progress_contract':'progress.json: phase training/inference, training_seconds; checkpoint.json: training_seconds,inference_seconds,blobs filename:sha256'}
    directory.mkdir(parents=True,exist_ok=True)
    # A retry must not accidentally accept artifacts left by a failed subprocess.
    for name in ('prediction.npy','case_ids.json','checkpoint.json','progress.json'):
        if (directory/name).exists(): (directory/name).unlink()
    write_json(directory/'request.json',request)
    started=time.monotonic()
    write_json(directory/'run.json',{'status':'running','identity':identity})
    invocation=command or [sys.executable,str(Path(__file__).with_name('baselines.py'))]
    environment=os.environ.copy()
    for name in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        environment[name]='2'
    try:
        with (directory/'adapter.log').open('w') as log:
            subprocess.run([*invocation,'--request',str((directory/'request.json').resolve())],
                           check=True,stdout=log,stderr=subprocess.STDOUT,env=environment)
        validate_prediction(directory,rows(files['test_cases']),request['n_genes'])
        checkpoint=read_json(directory/'checkpoint.json')
        checkpoint_artifacts(directory,checkpoint,require_blobs=command is not None)
        if checkpoint.get('identity') != identity:
            raise ValueError('checkpoint identity mismatch')
        run={'status':'completed','identity':identity,'elapsed_seconds':time.monotonic()-started,
             'training_seconds':checkpoint.get('training_seconds'), 'inference_seconds':checkpoint.get('inference_seconds'),
             'deterministic_reference':model in ('no_change','delta_mean'),
             'artifacts':{name:sha256(directory/name) for name in ('prediction.npy','case_ids.json','checkpoint.json','request.json')}}
        write_json(directory/'run.json',run)
        return run
    except Exception as exc:
        write_json(directory/'run.json',{'status':'failed','identity':identity,
                   'elapsed_seconds':time.monotonic()-started,'reason':str(exc)})
        raise


def freeze_predictions(input_root, output_root, models, commands=None):
    """Seal the exact 9-per-entry matrix only once every run is terminal."""
    commands=commands or {}
    runs=[]
    for model in models:
        for fold in range(3):
            for seed in SEEDS:
                directory=run_directory(output_root,model,fold,seed)
                try: run=read_json(directory/'run.json')
                except OSError as exc: raise ValueError('prediction barrier: missing run') from exc
                if run.get('status') not in TERMINAL:
                    raise ValueError('prediction barrier: nonterminal run')
                expected=run_identity(input_root,model,fold,seed,commands.get(model))
                if run.get('identity')!=expected:
                    raise ValueError('prediction barrier: stale identity')
                if run['status']=='completed': validate_run(input_root,output_root,model,fold,seed,commands.get(model))
                runs.append({'model':model,'fold':fold,'seed':seed,'status':run['status'],
                             'command':commands.get(model), 'run_sha256':sha256(directory/'run.json')})
    manifest={'input_sha256':sha256(Path(input_root)/'manifest.json'),'models':models,'runs':runs,
              'seeds':list(SEEDS),'folds':[0,1,2],'sealed':True}
    write_json(Path(output_root)/'predictions.json',manifest)
    return manifest


def verify_barrier(input_root, output_root):
    verify_input_files(input_root)
    manifest=read_json(Path(output_root)/'predictions.json')
    if not manifest.get('sealed') or manifest['input_sha256']!=sha256(Path(input_root)/'manifest.json'):
        raise ValueError('prediction barrier: input identity changed')
    expected={(m,f,s) for m in manifest['models'] for f in range(3) for s in SEEDS}
    actual=[(r['model'],r['fold'],r['seed']) for r in manifest['runs']]
    if set(actual)!=expected or len(actual)!=len(expected):
        raise ValueError('prediction barrier: incomplete or duplicate matrix')
    for entry in manifest['runs']:
        directory=run_directory(output_root,entry['model'],entry['fold'],entry['seed'])
        if entry['run_sha256']!=sha256(directory/'run.json') or entry['status'] not in TERMINAL:
            raise ValueError('prediction barrier: changed/nonterminal run')
        if entry['status']=='completed':
            validate_run(input_root,output_root,entry['model'],entry['fold'],entry['seed'],entry.get('command'))
    return manifest


def finite_mean(values, axis=None):
    values=np.asarray(values,dtype=float)
    count=np.sum(np.isfinite(values),axis=axis)
    return np.divide(np.nansum(values,axis=axis),count,
                     out=np.full(np.shape(count),np.nan,dtype=float),where=count>0)


def evaluate(input_root, output_root, workers=2):
    if workers not in (1,2): raise ValueError('evaluation workers must be 1 or 2')
    root=Path(input_root); output=Path(output_root)
    manifest=verify_barrier(root,output)
    # Fold truth/labels are read-only mmap and reused by all runs and both workers.
    truth={}
    for fold in range(3):
        directory=root/f'fold_{fold}'/'evaluator'
        truth[fold]=(rows(directory/'groups.csv'),
                     *(np.load(directory/name,mmap_mode='r') for name in ('truth.npy','positive.npy','negative.npy')))
    def score(entry):
        if entry['status']!='completed': return entry,[],None
        model,fold,seed=entry['model'],entry['fold'],entry['seed']
        directory=run_directory(output,model,fold,seed)
        prediction=np.load(directory/'prediction.npy',mmap_mode='r')
        cases=rows(root/f'fold_{fold}'/'test_cases.csv')
        groups,observed,positive,negative=truth[fold]
        group_rows=[]; curves=[]
        for index,group in enumerate(groups):
            selected=[i for i,case in enumerate(cases) if case['group_id']==group['group_id']]
            if not selected or len(selected)!=int(group['n_wells']): raise ValueError('test group membership mismatch')
            predicted=np.mean(prediction[selected],axis=0,dtype=np.float64)
            values,reasons,curve=score_group(observed[index],predicted,(positive[index],negative[index]))
            group_rows.append({'model':model,'fold':fold,'seed':seed,'group_id':group['group_id'],
                               **values, **{f'{key}_reason':reasons.get(key,'') for key in METRICS}})
            curves.append(curve)
        if sum(int(g['n_wells']) for g in groups)!=len(cases): raise ValueError('extra test cases')
        return entry,group_rows,finite_mean(curves,axis=0)
    with ThreadPoolExecutor(max_workers=workers) as pool: scored=list(pool.map(score,manifest['runs']))
    group_rows=[]; run_rows=[]; source_curves=[]; curve_ids=[]
    for entry,groups,curve in scored:
        group_rows.extend(groups)
        row={key:entry[key] for key in ('model','fold','seed','status')}
        row['deterministic_reference']=entry['model'] in ('no_change','delta_mean')
        for metric in METRICS:
            values=[g[metric] for g in groups]
            row[metric]=float(finite_mean(values))
            row[f'{metric}_valid_groups']=int(np.isfinite(values).sum())
        run_rows.append(row)
        if curve is not None:
            curve_ids.append(f"{entry['model']}:{entry['fold']}:{entry['seed']}")
            source_curves.append(curve)
    summary=[]; selections=[]; bands=[]
    for model in manifest['models']:
        model_runs=[r for r in run_rows if r['model']==model]
        completed=sum(r['status']=='completed' for r in model_runs)
        row={'model':model,'completed_runs':completed,'expected_runs':9,'complete_ranking':completed==9,
             'deterministic_reference':model in ('no_change','delta_mean')}
        for metric in METRICS:
            valid=np.array([r[metric] for r in model_runs]); valid=valid[np.isfinite(valid)]
            row[f'{metric}_mean']=float(valid.mean()) if len(valid) else np.nan
            row[f'{metric}_sample_sd']=float(valid.std(ddof=1)) if len(valid)>1 else np.nan
            row[f'{metric}_valid_runs']=len(valid)
        summary.append(row)
        eligible=[r for r in model_runs if np.isfinite(r[METRICS[0]]) and np.isfinite(r[METRICS[1]])]
        best=max(eligible,key=lambda r:(r[METRICS[0]]+r[METRICS[1]])/2) if eligible else None
        selections.append({'model':model,'fold':best['fold'] if best else '', 'seed':best['seed'] if best else '',
                           'selection_score':(best[METRICS[0]]+best[METRICS[1]])/2 if best else np.nan,
                           'label':'post-hoc best test; fixed fold/seed tie order'})
        model_curves=[c for identifier,c in zip(curve_ids,source_curves) if identifier.startswith(model+':')]
        for direction in range(2):
            values=np.array([c[direction] for c in model_curves])
            for index,recall in enumerate(RECALL_GRID):
                finite=values[:,index][np.isfinite(values[:,index])] if len(values) else []
                bands.append({'model':model,'direction':('positive','negative')[direction],'recall':recall,
                              'mean':float(np.mean(finite)) if len(finite) else np.nan,
                              'min':float(np.min(finite)) if len(finite) else np.nan,
                              'max':float(np.max(finite)) if len(finite) else np.nan,'valid_runs':len(finite)})
    write_csv(output/'per_group.csv',group_rows,fields=list(group_rows[0]) if group_rows else ['model','fold','seed','group_id',*METRICS])
    write_csv(output/'per_run.csv',run_rows)
    write_csv(output/'summary.csv',summary)
    write_csv(output/'best_runs.csv',selections)
    write_csv(output/'pr_bands.csv',bands)
    np.savez(output/'pr_sources.npz',recall=RECALL_GRID,identifiers=np.array(curve_ids),
             curves=np.array(source_curves) if source_curves else np.empty((0,2,len(RECALL_GRID))))
    paths=['per_group.csv','per_run.csv','summary.csv','best_runs.csv','pr_bands.csv','pr_sources.npz','predictions.json']
    write_json(output/'evaluation.json',{'artifacts':{p:sha256(output/p) for p in paths},
               'curve_convention':CURVE_CONVENTION,'workers':workers,'truth_cache':'once per fold, shared immutable mmap'})
    replay(output)
    return summary


def replay(output_root):
    output=Path(output_root)
    evidence=read_json(output/'evaluation.json')
    for name,digest in evidence['artifacts'].items():
        if sha256(output/name)!=digest: raise ValueError(f'replay checksum mismatch: {name}')
    from plots import render
    render(output)


def main():
    if len(sys.argv) > 1 and sys.argv[1] in ('evaluate', 'score', 'replay'):
        sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())))
        from stage2.evaluation.evaluate import main as shared_main
        if sys.argv[1] == 'evaluate':
            sys.argv[1] = 'score'
        sys.argv.extend(['--step', '3'])
        return shared_main()
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['synthetic','predict','seal','evaluate','replay','verify'])
    parser.add_argument('--input',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--model',default='no_change')
    parser.add_argument('--models',nargs='+',default=['no_change','delta_mean'])
    parser.add_argument('--fold',type=int,default=0); parser.add_argument('--seed',type=int,default=17)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--commands',type=Path,help='JSON mapping model names to isolated adapter command lists')
    args=parser.parse_args()
    commands=read_json(args.commands) if args.commands else {}
    if args.action=='replay': replay(args.output); return
    if args.input is None: parser.error('--input is required')
    if args.action=='synthetic':
        if not read_json(args.input/'manifest.json').get('synthetic'): raise ValueError('synthetic CLI requires marked fixture')
        for model in args.models:
            for fold in range(3):
                for seed in SEEDS: run_one(args.input,args.output,model,fold,seed,commands.get(model))
        freeze_predictions(args.input,args.output,args.models,commands)
        evaluate(args.input,args.output,args.workers)
    elif args.action=='predict': run_one(args.input,args.output,args.model,args.fold,args.seed,commands.get(args.model))
    elif args.action=='seal': freeze_predictions(args.input,args.output,args.models,commands)
    elif args.action=='evaluate': evaluate(args.input,args.output,args.workers)
    elif args.action=='verify': verify_barrier(args.input,args.output)

if __name__=='__main__': main()
