"""Full-matrix candidate acceptance and rollback-safe dose-model publication.

This module deliberately imports no model/matrix libraries at import time.
The expensive existing matrix acceptance and downstream replay run only after
all147 formal origins finish, under the repository cap and monitor.
"""
import json
import hashlib
import os
import shutil
import subprocess
import uuid
import csv
import tempfile
import sys
import argparse
import secrets
import math
from pathlib import Path

ACCEPTANCE='acceptance-cold-drug-held-dose-final.json'
STEP=Path(__file__).resolve().parents[1]
SEEDS=[17,42,2026]
ROLE='internal-validation-selected restored checkpoints; no new independent test'


def repository_root(start):
    """Locate a checkout, including worktrees whose .git marker is a file."""
    start=Path(start).resolve()
    for path in (start,*start.parents):
        if (path/'.git').exists():return path
    raise ValueError('repository root not found: no .git marker above '+str(start))


def read(path):return json.loads(Path(path).read_text())


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')


def digest(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):result.update(block)
    return result.hexdigest()


def file_snapshot(root,paths):
    root=Path(root)
    return {str(path):digest(root/path) for path in sorted(set(paths))}


def verify_file_snapshot(root,snapshot):
    if file_snapshot(root,snapshot)!=snapshot:raise ValueError('protected product preservation mismatch')


def capture_preservation(root):
    """Hash protected tracked bytes, acceptance manifests and native code only.

    This intentionally avoids rereading52GB of immutable data during training.
    Payload manifests pin their existing digests; the final offline health gate
    verifies large assets after the long trainer has stopped.
    """
    root=Path(root).resolve()
    protected=['stage2/01_data_raw','stage2/02_gene_axis','stage2/03_view','stage2/util',
               'stage2/04_model_view','stage3','pert_model/TranSiGen']
    raw=subprocess.check_output(['git','-C',str(root),'ls-files','-z','--',*protected])
    tracked=[p for p in raw.decode().split('\0') if p]
    manifests=[]
    for name in protected[:3]:
        manifests.extend(str(p.relative_to(root)) for p in (root/name/'output').rglob('manifest.json'))
    source=root/'pert_model/TranSiGen/source'
    native=[str(p.relative_to(root)) for p in source.rglob('*.py')]
    return dict(schema='tsg-dose-preservation-v1',scope='protected tracked bytes, accepted output manifests and immutable native Python source; large payload verification deferred to final offline health gate',files=file_snapshot(root,tracked+manifests+native),n_tracked=len(tracked),n_acceptance_manifests=len(manifests),n_native_source_files=len(native))


def verify_matrix(matrix):
    matrix=Path(matrix);accepted=read(matrix/ACCEPTANCE)
    if accepted.get('n_formal_fits')!=147:raise ValueError('publication requires complete147-fit acceptance')
    scenarios=['cold-drug','held-dose','final']
    methods=['no-dose','scaling','input-numeric','latent-numeric','input-embedding','latent-embedding','film']
    expected=[dict(scenario=s,fold=f,method=m,seed=t) for s in scenarios for f in range(1 if s=='final' else 3) for m in methods for t in [17,42,2026]]
    key=lambda value:json.dumps(value,sort_keys=True)
    if accepted.get('schema')!='tsg-dose-matrix-acceptance-v1' or accepted.get('scenarios')!=scenarios or sorted(map(key,accepted.get('runs',[])))!=sorted(map(key,expected)):
        raise ValueError('formal matrix unique coverage mismatch')
    if accepted.get('plan_sha256')!=digest(matrix/'plan.json') or read(matrix/'plan.json').get('runs')!=expected:raise ValueError('formal plan coverage/digest mismatch')
    for run in expected:
        path=matrix/'runs'/run['scenario']/run['method']/str(run['fold'])/str(run['seed'])
        terminal=read(path/'terminal.json')
        if terminal.get('status')!=('validation-complete' if run['scenario']=='final' else 'complete'):raise ValueError('incomplete or failed formal origin')
        manifest=read(path/'manifest.json')
        actual={p.name for p in path.iterdir() if p.is_file() and p.name!='manifest.json'}
        if manifest.get('schema')!='tsg-dose-fit-v1' or set(manifest.get('outputs',{}))!=actual:raise ValueError('formal origin manifest coverage mismatch')
        for name,value in manifest['outputs'].items():
            if Path(name).name!=name or digest(path/name)!=value:raise ValueError('formal origin payload digest mismatch')
    for terminal in (matrix/'attempts').rglob('terminal.json'):
        if read(terminal).get('status')=='failed':raise ValueError('failed scientific attempt blocks publication')
    return accepted


def assemble(matrix,inputs,candidate):
    """Copy only the final validation winner's original three checkpoints.

    Existing scientific acceptance must precede this structural packaging gate.
    No checkpoint deserialization, training or transform fitting occurs here.
    """
    matrix=Path(matrix);inputs=Path(inputs);candidate=Path(candidate)
    if candidate.exists():raise ValueError('release candidate must be an empty new path')
    verify_matrix(matrix)
    selection=matrix/'selections/final/0/receipt.json';receipt=read(selection)
    input_digest=digest(inputs/'manifest.json')
    if receipt.get('input_manifest_sha256')!=input_digest:raise ValueError('selection input identity mismatch')
    method=receipt['method']
    if method not in ['no-dose','scaling','input-numeric','latent-numeric','input-embedding','latent-embedding','film']:raise ValueError('unknown final selected method')
    candidate.mkdir(parents=True);(candidate/'checkpoints').mkdir()
    shutil.copy2(inputs/'genes.csv',candidate/'genes.csv')
    write(candidate/'input_contract.json',dict(read(inputs/'config.json'),input_manifest_sha256=input_digest))
    shutil.copy2(selection,candidate/'selection_receipt.json')
    counts={};environments={};origins={}
    for seed in SEEDS:
        origin=matrix/'runs/final'/method/'0'/str(seed)
        shutil.copy2(origin/'checkpoint.pt',candidate/'checkpoints'/f'{seed}.pt')
        counts[str(seed)]=read(origin/'terminal.json')['parameter_count']
        environments[str(seed)]=read(origin/'environment.json')
        origins[str(seed)]=dict(path=str(origin.resolve()),manifest_sha256=digest(origin/'manifest.json'))
    manifest=dict(schema='tsg-dose-release-v1',n_formal_fits=147,scenario_counts={'cold-drug':63,'held-dose':63,'final':21},method=method,seeds=SEEDS,role=ROLE,input_manifest_sha256=input_digest,parameter_counts=counts,environments=environments,origins=origins,matrix_acceptance_sha256=digest(matrix/ACCEPTANCE),matrix_plan_sha256=digest(matrix/'plan.json'),selection_sha256=digest(selection),outputs={str(p.relative_to(candidate)):digest(p) for p in sorted(candidate.rglob('*')) if p.is_file()})
    write(candidate/'manifest.json',manifest)
    return manifest


def csv_rows(path):
    with Path(path).open(newline='') as stream:return list(csv.DictReader(stream))


def write_csv(path,rows):
    if not rows:raise ValueError('empty plotted table: '+str(path))
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def held_dose_uncertainty(rows,repetitions=10000,seed=2262026):
    """Descriptive parent-block intervals from saved metrics, without reselection.

    Average the three seed metrics within each parent first. Every bootstrap
    draw carries that parent's complete dose/method/seed block. Undefined
    metrics stay excluded, with explicit parent and undefined-draw accounting.
    """
    import numpy as np
    if repetitions<1:raise ValueError('positive bootstrap count required')
    groups={}
    for row in rows:
        key=(row['method'],float(row['dose_um']));parent=str(row['drug_parent']);value=int(row['seed'])
        seeds=groups.setdefault(key,{}).setdefault(parent,{})
        if value in seeds:raise ValueError('duplicate held-dose seed metrics')
        seeds[value]=row
    if not groups:raise ValueError('empty held-dose parent metrics')
    parents=sorted(next(iter(groups.values())))
    arrays={};metrics=['pcc','spearman','rmse','r2']
    for key,blocks in sorted(groups.items()):
        if sorted(blocks)!=parents:raise ValueError('held-dose parent populations differ')
        for metric in metrics:
            values=[]
            for parent in parents:
                if set(blocks[parent])!=set(SEEDS):raise ValueError('incomplete held-dose seed metrics')
                scores=np.array([float(blocks[parent][s][metric]) if blocks[parent][s][metric] not in [None,''] else np.nan for s in SEEDS])
                if np.isinf(scores).any() or (np.isnan(scores).any() and not np.isnan(scores).all()):raise ValueError('inconsistent undefined seed metric')
                values.append(float(scores.mean()))
            arrays[key+(metric,)]=np.array(values)
    draws={key:[] for key in arrays};rng=np.random.default_rng(seed)
    for start in range(0,repetitions,128):
        index=rng.integers(len(parents),size=(min(128,repetitions-start),len(parents)))
        for key,values in arrays.items():
            samples=values[index];defined=np.isfinite(samples);counts=defined.sum(axis=1)
            totals=np.where(defined,samples,0).sum(axis=1)
            means=np.full(len(index),np.nan);np.divide(totals,counts,out=means,where=counts>0)
            draws[key].append(means)
    result=[]
    for (method,dose,metric),values in arrays.items():
        sampled=np.concatenate(draws[(method,dose,metric)]);defined=values[np.isfinite(values)]
        valid=sampled[np.isfinite(sampled)]
        low,high=np.quantile(valid,[.025,.975]) if len(valid) else (None,None)
        result.append(dict(method=method,dose_um=dose,metric=metric,estimate=float(defined.mean()) if len(defined) else None,low=float(low) if low is not None else None,high=float(high) if high is not None else None,n_parents=len(parents),n_defined_parents=len(defined),n_excluded_undefined=len(parents)-len(defined),n_seeds=3,repetitions=repetitions,statistical_seed=seed,n_bootstrap_undefined=repetitions-len(valid),inference='descriptive 95% percentile parent-block interval; all doses and seeds carried together; no method reselection',role='auxiliary descriptive; three measured doses only'))
    return result


def export_metric_grid(reports,destination):
    """Record the same four indicator panels as step02, from saved metrics only."""
    result=[]
    for scenario in ['cold-drug','held-dose']:
        source=reports/scenario
        parents=csv_rows(source/'per_parent.csv')
        for summary in csv_rows(source/'summary.csv'):
            records=[r for r in parents if r['method']==summary['method']]
            for metric in ['pcc','spearman','rmse','r2']:
                defined=[r for r in records if r[metric] not in ['',None] and math.isfinite(float(r[metric]))]
                result.append(dict(scenario=scenario,method=summary['method'],metric=metric,value=summary[metric],
                    n_parents=len(records),n_defined=len(defined),n_undefined=len(records)-len(defined),
                    n_constant_truth=sum(int(r.get('truth_constant_count',0)) for r in records),
                    n_constant_prediction=sum(int(r.get('prediction_constant_count',0)) for r in records),
                    n_undefined_parent_dose=sum(int(r.get('undefined_r2_count' if metric=='r2' else 'undefined_pcc_count',0)) for r in records) if metric!='rmse' else 0,
                    zero=0,direction='lower is better' if metric=='rmse' else 'higher is better',
                    role='primary selected procedure' if summary['method']=='selected' else ('auxiliary fixed comparison' if scenario=='held-dose' else 'descriptive fixed candidate'),
                    source=str(Path('reports')/scenario/'summary.csv'),source_sha256=digest(source/'summary.csv'),
                    counts_source_sha256=digest(source/'per_parent.csv')))
    write_csv(destination/'metrics_values.csv',result)


def export_plot_style(destination):
    """Snapshot the pure shared Stage2 palette so isolated rendering reads CSV only."""
    import importlib.util
    path=STEP.parent/'util/output/api/stage2_util/palette.py'
    spec=importlib.util.spec_from_file_location('dose_report_palette',path)
    palette=importlib.util.module_from_spec(spec);spec.loader.exec_module(palette)
    hues=dict(zip(['no-dose','scaling','input-numeric','latent-numeric','input-embedding','latent-embedding','film','selected'],['gray','orange','blue','green','purple','magenta','blue','gray']))
    values=[dict(role='method',name=method,line=palette.color(hue,0 if method in ['film','selected'] else 2),fill=palette.color(hue,4 if method in ['film','selected'] else 6),source_sha256=digest(path)) for method,hue in hues.items()]
    for index in range(9):
        values.append(dict(role='run',name=index,train=palette.color('blue',palette.run_shade(index)),validation=palette.color('green',palette.run_shade(index)),source_sha256=digest(path)))
    write_csv(destination/'plot_style.csv',values)


def export_plot_data(matrix,destination):
    """Keep saved points and add descriptive uncertainty, never rescore predictions."""
    matrix=Path(matrix);destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    reports=matrix/'reports';selections=matrix/'selections'
    export_metric_grid(reports,destination);export_plot_style(destination)
    final=read(selections/'final/0/receipt.json')['method']
    references={f'primary_fold_{f}_selected_method':read(selections/'cold-drug'/str(f)/'receipt.json')['method'] for f in range(3)}
    references['final_selected_method']=final
    for source,target in [('summary.csv','primary_values.csv'),('primary_paired.csv','primary_paired.csv'),('descriptive_paired.csv','primary_descriptive.csv'),('per_dose.csv','primary_dose_values.csv')]:
        shutil.copy2(reports/'cold-drug'/source,destination/target)
    intervals=held_dose_uncertainty(csv_rows(reports/'held-dose/per_parent_dose_seed.csv'))
    intervals=[dict(row,**references,evaluation={.2:'lower extrapolation',1.:'interpolation',10.:'upper extrapolation'}[row['dose_um']]) for row in intervals]
    write_csv(destination/'held_dose_uncertainty.csv',intervals)
    lookup={(row['method'],row['dose_um'],row['metric']):row for row in intervals}
    held=[]
    for row in csv_rows(reports/'held-dose/per_dose.csv'):
        row=dict(row,**references,evaluation={.2:'lower extrapolation',1.:'interpolation',10.:'upper extrapolation'}[float(row['dose_um'])])
        for metric in ['pcc','spearman','rmse','r2']:
            interval=lookup[(row['method'],float(row['dose_um']),metric)]
            if metric in row and interval['estimate'] is not None and not math.isclose(float(row[metric]),interval['estimate'],rel_tol=1e-12,abs_tol=1e-12):raise ValueError('held-dose point estimate differs from saved parent metrics')
            for key in ['low','high','n_defined_parents','n_excluded_undefined','n_bootstrap_undefined']:row[metric+'_'+key]=interval[key]
        row.update(repetitions=10000,statistical_seed=2262026,inference=interval['inference'])
        held.append(row)
    write_csv(destination/'held_dose_values.csv',held)
    curves=[];resources=[];validations=[]
    for scenario in ['cold-drug','held-dose','final']:
        curves.extend(csv_rows(reports/scenario/'learning_curves.csv'))
        resources.extend(csv_rows(reports/scenario/'fit_resources.csv'))
        for fold in range(1 if scenario=='final' else 3):
            receipt=read(selections/scenario/str(fold)/'receipt.json')
            for row in csv_rows(selections/scenario/str(fold)/'validation_table.csv'):
                validations.append(dict(row,scenario=scenario,fold=fold,selected_method=receipt['method'],role='fixed comparison only' if scenario=='held-dose' else 'validation selection'))
    write_csv(destination/'learning_curves_all.csv',curves)
    write_csv(destination/'resources_all.csv',resources)
    write_csv(destination/'final_validation_values.csv',validations)
    return {p.name:digest(p) for p in sorted(destination.glob('*.csv'))}


def render_isolated(data,destination,python):
    """Give the renderer only its script and CSV bytes, then retain the complete metric and learning-curve PNG set."""
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='tsg-dose-csv-only-') as tmp:
        isolated=Path(tmp);(isolated/'data').mkdir()
        shutil.copy2(STEP/'figure/script/render.py',isolated/'render.py')
        for path in Path(data).glob('*.csv'):shutil.copy2(path,isolated/'data'/path.name)
        subprocess.run([str(python),'-I',str(isolated/'render.py'),'--data',str(isolated/'data'),'--output',str(isolated/'figure')],cwd=isolated,env=dict(os.environ,MPLCONFIGDIR=str(isolated/'mpl-cache')),check=True)
        expected={'primary.png','primary_dose.png','held_dose.png','final_validation.png','learning_curves.png','resources.png','metrics.png','paired_differences.png','learning_curve.png','learning_curve_cold-drug.png','learning_curve_held-dose.png','learning_curve_final.png'}
        if {p.name for p in (isolated/'figure').glob('*.png')}!=expected:raise ValueError('incomplete CSV-only figure reproduction')
        for name in expected:shutil.copy2(isolated/'figure'/name,destination/name)
    return {name:digest(destination/name) for name in sorted(expected)}


def verify_downstream(candidate,inputs,matrix,source,receipt_path):
    """Replay all validation rows through the standalone reader, without labels.

    Run only in a fresh, capped process after formal fitting has finished. Exact
    full-row order preserves sampled native batch64 randomness at every seed.
    """
    sys.path.insert(0,str(STEP/'output/api'))
    import numpy as np
    from tsg_dose.inputs import Inputs
    from tsg_dose.release import Release
    candidate=Path(candidate);matrix=Path(matrix);data=Inputs(inputs)
    release=Release(candidate,source);rows=data.rows('final',0,'validation')
    ids=data.obs.index[rows].astype(str).tolist();gene_ids=data.genes.gene_index.astype(str).tolist();proofs=[]
    for seed in SEEDS:
        origin=matrix/'runs/final'/release.manifest['method']/'0'/str(seed)
        config=read(origin/'config.json')
        if config['validation_order']!=rows.tolist():raise ValueError('downstream validation row order mismatch')
        prediction=release.predict(seed,controls=data.control[rows],morgan=data.fingerprints(rows),doses=data.doses(rows),gene_ids=gene_ids,profile_ids=ids,sampling_seed=config['validation_seed'])
        expected=np.load(origin/'validation_prediction.npy',mmap_mode='r')
        if prediction.index.astype(str).tolist()!=ids or prediction.columns.tolist()!=gene_ids or prediction.to_numpy().dtype!=np.float32 or not np.array_equal(prediction.to_numpy(),expected):raise ValueError('standalone downstream exact replay mismatch')
        proofs.append(dict(seed=seed,n_rows=len(rows),n_genes=len(gene_ids),prediction_sha256=digest(origin/'validation_prediction.npy'),checkpoint_sha256=digest(candidate/'checkpoints'/f'{seed}.pt')))
    result=dict(schema='tsg-dose-downstream-replay-v1',exact=True,pid=os.getpid(),process_token=secrets.token_hex(32),candidate_manifest_sha256=digest(candidate/'manifest.json'),proofs=proofs)
    write(receipt_path,result);return result


def retain_audit(matrix,log):
    """Keep compact origin and choice evidence beside the acceptance receipt."""
    matrix=Path(matrix);log=Path(log);origins=[];choices=[];attempts=[]
    for run in read(matrix/ACCEPTANCE)['runs']:
        path=matrix/'runs'/run['scenario']/run['method']/str(run['fold'])/str(run['seed'])
        terminal=read(path/'terminal.json');config=read(path/'config.json')
        origins.append(dict(run,**{k:terminal[k] for k in ['status','best_epoch','epochs','parameter_count','elapsed_seconds','peak_rss_kib']},transform=json.dumps(config['transform'],sort_keys=True),checkpoint_sha256=digest(path/'checkpoint.pt'),config_sha256=digest(path/'config.json'),environment_sha256=digest(path/'environment.json'),manifest_sha256=digest(path/'manifest.json'),validation_replay=json.dumps(terminal['validation_replay'],sort_keys=True),test_replay=json.dumps(terminal.get('test_replay'),sort_keys=True)))
    for path in sorted((matrix/'selections').glob('*/*/receipt.json')):choices.append(dict(path=str(path.relative_to(matrix)),sha256=digest(path),receipt=read(path)))
    for path in sorted((matrix/'attempts').rglob('terminal.json')):attempts.append(dict(path=str(path.relative_to(matrix)),sha256=digest(path),terminal=read(path)))
    write_csv(log.with_name(log.stem+'-origins.csv'),origins)
    write(log.with_name(log.stem+'-choices.json'),dict(choices=choices,retained_attempts=attempts))


def publish(args):
    """Complete scientific acceptance, replay and figures before stable install."""
    # Existing scientific verifier covers registrations, configuration, all
    # predictions/metrics and independent replay evidence. No training allowed.
    subprocess.run([args.python,str(STEP/'execution/run_matrix.py'),'accept','--inputs',str(args.inputs),'--candidate',str(args.matrix),'--source',str(args.source),'--device',args.device],check=True)
    assemble(args.matrix,args.inputs,args.candidate)
    downstream=args.log.with_name(args.log.stem+'-downstream.json')
    if downstream.exists():raise ValueError('downstream receipt path already exists')
    subprocess.run([args.python,str(Path(__file__).resolve()),'downstream','--inputs',str(args.inputs),'--matrix',str(args.matrix),'--candidate',str(args.candidate),'--source',str(args.source),'--log',str(downstream)],check=True)
    proof=read(downstream)
    if proof.get('exact') is not True or proof.get('pid')==os.getpid() or not proof.get('process_token') or proof.get('candidate_manifest_sha256')!=digest(args.candidate/'manifest.json') or [r['seed'] for r in proof.get('proofs',[])]!=SEEDS:raise ValueError('missing fresh-process downstream acceptance')
    args.figures.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.tsg-dose-figures-',dir=args.figures.parent) as tmp:
        staged=Path(tmp)/'figure'
        if args.figures.exists():shutil.copytree(args.figures,staged)
        data_hashes=export_plot_data(args.matrix,staged/'data')
        figure_hashes=render_isolated(staged/'data',staged,args.python)
        before=read(args.preservation);after=capture_preservation(args.root)
        if before['files']!=after['files']:raise ValueError('protected product preservation mismatch')
        retain_audit(args.matrix,args.log)
        receipt=dict(schema='tsg-dose-publication-acceptance-v1',status='passed',n_formal_fits=147,downstream_exact_replay=True,csv_only_figures=True,protected_products_unchanged=True,candidate_manifest_sha256=digest(args.candidate/'manifest.json'),matrix_acceptance_sha256=digest(args.matrix/ACCEPTANCE),downstream_receipt_sha256=digest(downstream),figure_csv_sha256=data_hashes,figure_png_sha256=figure_hashes,preservation=dict(n_files=len(before['files']),before_sha256=digest(args.preservation)),role=ROLE)
        receipt['installation']=install(args.candidate,args.destination,receipt,figures_candidate=staged,figures_destination=args.figures)
    write(args.log,receipt);return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['publish','downstream'])
    for name in ['inputs','matrix','candidate','source','log']:parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--destination',type=Path);parser.add_argument('--figures',type=Path,default=STEP/'figure')
    parser.add_argument('--preservation',type=Path,default=STEP/'log/239-preservation-before.json')
    parser.add_argument('--root',type=Path);parser.add_argument('--python',default=sys.executable);parser.add_argument('--device',default='cuda')
    args=parser.parse_args()
    if args.root is None:args.root=repository_root(STEP)
    for name in ['inputs','matrix','candidate','source','log','figures','preservation','root']:setattr(args,name,getattr(args,name).resolve())
    if args.action=='downstream':verify_downstream(args.candidate,args.inputs,args.matrix,args.source,args.log)
    else:
        if args.destination is None:parser.error('publish requires --destination')
        args.destination=args.destination.resolve();publish(args)


def install(candidate,destination,acceptance,*,figures_candidate=None,figures_destination=None):
    """Install only an accepted candidate; rollback any failed replacement.

    The acceptance receipt is produced after actual downstream replay and
    isolated CSV-only figure rendering. Model API code is never replaced.
    """
    candidate=Path(candidate);destination=Path(destination)
    if acceptance.get('status')!='passed' or acceptance.get('n_formal_fits')!=147 or acceptance.get('downstream_exact_replay') is not True or acceptance.get('csv_only_figures') is not True or acceptance.get('candidate_manifest_sha256')!=digest(candidate/'manifest.json'):
        raise ValueError('candidate has no complete downstream/figure acceptance')
    if (figures_candidate is None)!=(figures_destination is None):raise ValueError('figure candidate and destination must be supplied together')
    pairs=[(candidate,destination)]
    if figures_candidate is not None:pairs.append((Path(figures_candidate),Path(figures_destination)))
    paths=[p.resolve() for pair in pairs for p in pair]
    if len(set(paths))!=len(paths) or any(p.is_symlink() for pair in pairs for p in pair):raise ValueError('unsafe release installation destination')
    changes=[]
    try:
        for source,target in pairs:
            target.parent.mkdir(parents=True,exist_ok=True)
            backup=target.with_name('.'+target.name+'.backup-'+uuid.uuid4().hex)
            change=dict(source=source,target=target,backup=backup,had_previous=target.exists(),installed=False)
            changes.append(change)
            if change['had_previous']:os.replace(target,backup)
            os.replace(source,target);change['installed']=True
        if digest(destination/'manifest.json')!=acceptance['candidate_manifest_sha256']:raise ValueError('installed release digest changed')
    except BaseException:
        for change in reversed(changes):
            if change['installed'] and change['target'].exists():os.replace(change['target'],change['source'])
            if change['had_previous'] and change['backup'].exists():os.replace(change['backup'],change['target'])
        raise
    for change in changes:
        if change['backup'].exists():shutil.rmtree(change['backup'])
    return dict(destination=str(destination.resolve()),manifest_sha256=acceptance['candidate_manifest_sha256'],replaced_previous=changes[0]['had_previous'],figures_destination=str(Path(figures_destination).resolve()) if figures_destination is not None else None)


if __name__=='__main__':main()
