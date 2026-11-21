"""Export and accept round1 reports from saved metrics; never fits or predicts."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import tempfile

STEP=Path(__file__).resolve().parents[1]
METHODS=['no-dose','scaling','input-numeric','latent-numeric','input-embedding','latent-embedding','film']
SCENARIOS=['cold-drug','held-dose']
SEEDS=[17,42,2026]
METRICS=['pcc','spearman','rmse','r2']
EXPECTED={'round1_metrics.png','round1_learning_curves.png','round1_metrics_cold_drug_dose.png'}

def digest(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):value.update(block)
    return value.hexdigest()

def rows(path):
    with Path(path).open(newline='') as stream:return list(csv.DictReader(stream))

def write_csv(path,values):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(values[0]));writer.writeheader();writer.writerows(values)

def number(value):
    return float(value) if value not in ['',None] else float('nan')

def finite_mean(values):
    finite=[v for v in values if math.isfinite(v)]
    return statistics.mean(finite) if finite else float('nan')

def export_run_metrics(matrix,destination):
    matrix=Path(matrix);destination=Path(destination);output=[]
    for scenario in SCENARIOS:
        for method_index,method in enumerate(METHODS):
            for fold in range(3):
                selection=matrix/f'selections/{scenario}/{fold}/receipt.json'
                selected=json.loads(selection.read_text())['method']
                if selected not in METHODS:raise ValueError('invalid sealed selected method')
                for seed_index,seed in enumerate(SEEDS):
                    source=Path(f'runs/{scenario}/{method}/{fold}/{seed}/test_per_group.csv')
                    values=rows(matrix/source)
                    if not values:raise ValueError('empty saved parent-dose metrics')
                    identities=[(r['drug_parent'],float(r['dose_um'])) for r in values]
                    if len(set(identities))!=len(identities):raise ValueError('duplicate saved parent-dose metric')
                    for metric in METRICS:
                        parents=defaultdict(list)
                        for row in values:parents[row['drug_parent']].append(number(row[metric]))
                        parent_means=[finite_mean(part) for part in parents.values()]
                        raw=[v for part in parents.values() for v in part]
                        point=dict(scenario=scenario,method=method,actual_method=method,fold=fold,seed=seed,metric=metric,value=finite_mean(parent_means),
                                   x=method_index+(fold*3+seed_index-4)*.045,
                                   n_parent_dose=len(raw),n_nonfinite_parent_dose=sum(not math.isfinite(v) for v in raw),
                                   n_parents=len(parents),n_undefined_parents=sum(not math.isfinite(v) for v in parent_means),
                                   n_constant_truth=sum(str(r['constant_truth']).lower() in ['true','1'] for r in values),n_constant_prediction=sum(str(r['constant_prediction']).lower() in ['true','1'] for r in values),
                                   source=str(source),source_sha256=digest(matrix/source),selection_source='',selection_sha256='')
                        output.append(point)
                        if scenario=='cold-drug' and method==selected:
                            output.append(dict(point,method='selected',x=7+(fold*3+seed_index-4)*.045,selection_source=str(selection.relative_to(matrix)),selection_sha256=digest(selection)))
    output.sort(key=lambda r:(SCENARIOS.index(r['scenario']),(METHODS+['selected']).index(r['method']),METRICS.index(r['metric']),r['fold'],SEEDS.index(r['seed'])))
    write_csv(destination/'round1_runs.csv',output)
    summaries=[]
    for scenario in SCENARIOS:
        for method in METHODS+(['selected'] if scenario=='cold-drug' else []):
            for metric in METRICS:
                part=[r for r in output if (r['scenario'],r['method'],r['metric'])==(scenario,method,metric)]
                finite=[r['value'] for r in part if math.isfinite(r['value'])]
                summaries.append(dict(scenario=scenario,method=method,x=(METHODS+['selected']).index(method),metric=metric,source='round1_runs.csv',source_sha256=digest(destination/'round1_runs.csv'),mean=finite_mean(finite),sd=statistics.stdev(finite) if len(finite)>1 else float('nan'),n_total=len(part),n_defined=len(finite),n_undefined=len(part)-len(finite),n_nonfinite_parent_dose=sum(r['n_nonfinite_parent_dose'] for r in part),n_undefined_parents=sum(r['n_undefined_parents'] for r in part),aggregation='equal run mean; sample SD ddof=1; finite values only; no zero substitution'))
    write_csv(destination/'round1_summary.csv',summaries)

def export_coordinates(matrix,destination):
    matrix=Path(matrix);destination=Path(destination)
    source=Path('reports/cold-drug/per_dose.csv');values=rows(matrix/source)
    expected={(m,d) for m in METHODS+['selected'] for d in [.2,1.,10.]}
    if {(r['method'],float(r['dose_um'])) for r in values}!=expected or len(values)!=24:
        raise ValueError('cold-drug dose coverage mismatch')
    output=[]
    for row in values:
        for metric in METRICS:
            output.append(dict(method=row['method'],dose_um=row['dose_um'],metric=metric,value=row[metric],x_log10=math.log10(float(row['dose_um'])),n_parents=row['n_parents'],n_seeds=row['n_seeds'],aggregation='parent-equal within seed; equal seed OOF',source=str(source),source_sha256=digest(matrix/source)))
    write_csv(destination/'round1_dose.csv',output)
    curves=[]
    expected_runs={(m,str(f),str(s)) for m in METHODS for f in range(3) for s in SEEDS}
    for scenario in SCENARIOS:
        source=Path(f'reports/{scenario}/learning_curves.csv');part=rows(matrix/source);source_digest=digest(matrix/source)
        if {(r['method'],r['fold'],r['seed']) for r in part}!=expected_runs:raise ValueError('curve run coverage mismatch')
        identities=[(r['method'],r['fold'],r['seed'],r['epoch']) for r in part]
        if len(identities)!=len(set(identities)):raise ValueError('duplicate curve coordinate')
        for row in part:
            if row['scenario']!=scenario:raise ValueError('curve scenario mismatch')
            if any(not math.isfinite(number(row[k])) for k in ['train_loss','validation_mse']) or number(row['validation_mse'])<=0:
                raise ValueError('nonfinite or nonpositive curve coordinate')
            curves.append(dict(row,source=str(source),source_sha256=source_digest))
    write_csv(destination/'round1_curves.csv',curves)

def historical_report_root():
    archive = STEP / 'execution/historical/report_126'
    return archive if archive.is_dir() else STEP / 'figure'


def render_isolated(data,destination,python):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='round1-csv-only-') as tmp:
        isolated=Path(tmp);shutil.copytree(data,isolated/'data');shutil.copy2(historical_report_root()/'script/render.py',isolated/'render.py')
        subprocess.run([str(python),'-I',str(isolated/'render.py'),'--data',str(isolated/'data'),'--output',str(isolated/'rendered')],check=True,env=dict(os.environ,MPLCONFIGDIR=str(isolated/'mpl')))
        rendered={p.name:digest(p) for p in (isolated/'rendered').glob('*.png')}
        if set(rendered)!=EXPECTED:raise ValueError('round1 figure filename coverage mismatch')
        for name in rendered:shutil.copy2(isolated/'rendered'/name,destination/name)
    return rendered

def snapshot(directory):
    directory=Path(directory)
    return {str(p.relative_to(directory)):digest(p) for p in sorted(directory.rglob('*')) if p.is_file()}

def install(candidate,destination,receipt):
    """Install only accepted report bytes; roll back a filesystem failure."""
    candidate=Path(candidate);destination=Path(destination)
    if (destination.resolve() == (STEP/'figure').resolve()
            and (STEP/'execution/historical/report_126').is_dir()):
        raise ValueError('canonical active figure now belongs to the dual-source release')
    if not (receipt.get('status')=='passed' and receipt.get('n_formal_fits')==126 and receipt.get('csv_only_figures') and receipt.get('visual_inspection')):
        raise ValueError('report acceptance incomplete')
    if snapshot(candidate)!=receipt.get('files') or {p.name for p in candidate.glob('*.png')}!=EXPECTED:
        raise ValueError('accepted candidate changed')
    backup=destination.with_name(destination.name+'.round1-backup')
    if backup.exists():raise ValueError('report backup already exists')
    existed=destination.exists()
    if existed:os.replace(destination,backup)
    try:os.replace(candidate,destination)
    except BaseException:
        if existed:os.replace(backup,destination)
        raise
    if existed:shutil.rmtree(backup)


def verify_scientific_receipt(matrix,receipt_path):
    matrix=Path(matrix);receipt=json.loads(Path(receipt_path).read_text())
    expected=[dict(scenario=scenario,fold=fold,method=method,seed=seed) for scenario in SCENARIOS for fold in range(3) for method in METHODS for seed in SEEDS]
    key=lambda r:(r['scenario'],r['fold'],r['method'],r['seed'])
    if receipt.get('schema')!='tsg-dose-round1-acceptance-v1' or receipt.get('scenarios')!=SCENARIOS or receipt.get('n_formal_fits')!=126 or sorted(map(key,receipt.get('runs',[])))!=sorted(map(key,expected)):
        raise ValueError('complete126-origin round1 scientific acceptance required')
    if receipt.get('plan_sha256')!=digest(matrix/'plan.json'):raise ValueError('accepted scientific plan changed')
    choices={f'selections/{scenario}/{fold}/receipt.json' for scenario in SCENARIOS for fold in range(3)}
    if set(receipt.get('selection_sha256',{}))!=choices:
        raise ValueError('accepted selection receipt coverage mismatch')
    for name,value in receipt['selection_sha256'].items():
        if digest(matrix/name)!=value:raise ValueError('accepted selection receipt changed')
    manifests={f"runs/{r['scenario']}/{r['method']}/{r['fold']}/{r['seed']}/manifest.json" for r in expected}
    if set(receipt.get('origin_manifest_sha256',{}))!=manifests:
        raise ValueError('accepted origin manifest coverage mismatch')
    for name,value in receipt['origin_manifest_sha256'].items():
        if digest(matrix/name)!=value:raise ValueError('accepted origin manifest changed')
    reports={f'reports/{scenario}/learning_curves.csv' for scenario in SCENARIOS}|{'reports/cold-drug/per_dose.csv'}
    if set(receipt.get('report_source_sha256',{}))!=reports:
        raise ValueError('accepted scientific report coverage mismatch')
    for name,value in receipt['report_source_sha256'].items():
        if digest(matrix/name)!=value:raise ValueError('accepted scientific report changed')
    # Bind the saved scalar report inputs to the origin manifests already replayed
    # by scientific acceptance, without generating any new predictions.
    for run in expected:
        origin=matrix/f"runs/{run['scenario']}/{run['method']}/{run['fold']}/{run['seed']}"
        manifest=json.loads((origin/'manifest.json').read_text())
        if manifest.get('outputs',{}).get('test_per_group.csv')!=digest(origin/'test_per_group.csv'):
            raise ValueError('accepted origin report metric changed')
    return receipt


def prepare(matrix,candidate,previous,scientific_receipt,python):
    matrix=Path(matrix);candidate=Path(candidate);previous=Path(previous)
    scientific=verify_scientific_receipt(matrix,scientific_receipt)
    if candidate.exists():raise ValueError('report candidate must be new')
    candidate.mkdir(parents=True)
    try:
        shutil.copytree(previous/'data',candidate/'data');shutil.copytree(STEP/'figure/script',candidate/'script')
        export_run_metrics(matrix,candidate/'data');export_coordinates(matrix,candidate/'data')
        rendered=render_isolated(candidate/'data',candidate,python)
        sources={r['source']:r['source_sha256'] for name in ['round1_runs.csv','round1_curves.csv','round1_dose.csv'] for r in rows(candidate/'data'/name)}
        sources.update(scientific['selection_sha256'])
        sources.update(scientific['origin_manifest_sha256'])
        sources.update(scientific['report_source_sha256'])
        for name,value in sources.items():
            if digest(matrix/name)!=value:raise ValueError('report source changed during export')
        old={p.name:digest(p) for p in previous.glob('*.png')}
        new_roles={'metrics.png':'round1_metrics.png','primary.png':'round1_metrics.png','primary_dose.png':'round1_metrics_cold_drug_dose.png','learning_curve.png':'round1_learning_curves.png','learning_curves.png':'round1_learning_curves.png','learning_curve_cold-drug.png':'round1_learning_curves.png','learning_curve_held-dose.png':'round1_learning_curves.png','held_dose.png':'round1_metrics.png'}
        return dict(schema='tsg-dose-round1-report-acceptance-v1',status='candidate-verified',n_formal_fits=126,csv_only_figures=True,visual_inspection=False,
                    scientific_receipt_sha256=digest(scientific_receipt),matrix=str(matrix.resolve()),sources=sources,
                    files=snapshot(candidate),figures=rendered,
                    superseded_figures={name:dict(sha256=value,role=new_roles.get(name,'historical scientific evidence only; no active PNG')) for name,value in old.items()},
                    undefined_policy='Finite doses averaged within parent, finite parents within run, finite run values in report mean/SD; all excluded counts explicit; no zero substitution',
                    aggregation='run-equal mean and sample SD ddof=1; dose summaries retain parent-equal/equal-seed OOF',training_or_reprediction=False)
    except BaseException:
        shutil.rmtree(candidate)
        raise


def accept(candidate,receipt):
    candidate=Path(candidate)
    if receipt.get('status')!='candidate-verified' or snapshot(candidate)!=receipt.get('files'):
        raise ValueError('candidate content changed before acceptance')
    for name,value in receipt['sources'].items():
        if digest(Path(receipt['matrix'])/name)!=value:raise ValueError('saved report source changed')
    return dict(receipt,status='passed',visual_inspection=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    stage=sub.add_parser('prepare')
    for name in ['matrix','candidate','previous','scientific-receipt','receipt']:stage.add_argument('--'+name,type=Path,required=True)
    stage.add_argument('--python',default=sys.executable)
    approve=sub.add_parser('accept');approve.add_argument('--candidate',type=Path,required=True);approve.add_argument('--receipt',type=Path,required=True);approve.add_argument('--visual-inspection',action='store_true',required=True)
    publish=sub.add_parser('install')
    for name in ['candidate','destination','receipt']:publish.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    if args.command=='prepare':
        receipt=prepare(args.matrix,args.candidate,args.previous,args.scientific_receipt,args.python)
    else:
        receipt=json.loads(args.receipt.read_text())
        if args.command=='accept':receipt=accept(args.candidate,receipt)
        else:install(args.candidate,args.destination,receipt);return
    args.receipt.parent.mkdir(parents=True,exist_ok=True);args.receipt.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':main()
