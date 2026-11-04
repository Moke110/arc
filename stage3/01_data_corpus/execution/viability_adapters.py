#!/usr/bin/env python3
"""Freeze reviewed figure curation or emit an offline source-evidence bundle."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
HERE=Path(__file__).resolve().parent
SOURCES={'linn_2024':{'item':'Linn-2024','curation':'linn_v1.json','image':'ijms-25-01427-g003.jpg','method':'MTT'}}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def prepare(source):
    """Explicit curation setup; ordinary bundle generation remains read-only."""
    config=SOURCES[source]
    item=ROOT/'pheno_data'/config['item']
    tracked=HERE/'curation'/config['curation']
    destination=item/'local'/config['curation']
    destination.parent.mkdir(exist_ok=True)
    if destination.exists() and destination.read_bytes()!=tracked.read_bytes():
        raise ValueError('Versioned curation differs; publish a new version')
    destination.write_bytes(tracked.read_bytes())
    path=item/'author_log.json'
    provenance=json.loads(path.read_text())
    entry={'path':'local/'+config['curation'],'sha256':sha(destination),'bytes':destination.stat().st_size,
           'origin':str(tracked.relative_to(ROOT))+'; reviewed author-summary figure digitization, version1',
           'source_image_sha256':sha(item/'source'/config['image'])}
    prior=next((r for r in provenance['local_additions'] if r['path']==entry['path']),None)
    if prior is not None and prior!=entry: raise ValueError('Declared curation provenance changed')
    if prior is None: provenance['local_additions'].append(entry)
    path.write_text(json.dumps(provenance,indent=2)+'\n')

def source_context(source):
    config=SOURCES[source]
    item=ROOT/'pheno_data'/config['item']
    provenance=json.loads((item/'author_log.json').read_text())
    assets=[]
    for entry in provenance['sources']:
        assets.append({'path':str((item/'source'/entry['path']).relative_to(ROOT)),'sha256':entry['sha256']})
    for entry in provenance['local_additions']:
        assets.append({'path':str((item/entry['path']).relative_to(ROOT)),'sha256':entry['sha256']})
    for entry in assets:
        if sha(ROOT/entry['path'])!=entry['sha256']: raise ValueError('Frozen asset mismatch '+entry['path'])
    curation_path=item/'local'/config['curation']
    if curation_path.read_bytes()!=(HERE/'curation'/config['curation']).read_bytes():
        raise ValueError('Curation snapshot differs from versioned extraction')
    data=json.loads(curation_path.read_text())
    context={'source':source,'item':'pheno_data/'+config['item'],'test_only':False,'status':'archived','assets':assets}
    return context,data,config

def condition(cell,drugtime,hypoxia,reoxygenation,**extra):
    return {'condition_cell_type':cell,'condition_drug_time':drugtime,'condition_hypoxia_time':hypoxia,
            'condition_reoxygen_time':reoxygenation,'condition_drug_schedule':'direct' if hypoxia==0 else 'pretreatment; persistence during injury unknown',**extra}

def linn_records(source,data,config):
    records=[]
    for index,point in enumerate(data['points'],1):
        identity=f'{source}-Fig3b-{index:02}'
        cond=condition(point['cell'],48,0,0,condition_species='Homo sapiens',condition_model='MUi019-A' if point['cell']=='hiPSC-CM' else 'AC16',condition_infection='absent in cytotoxicity experiment',condition_vehicle='not explicitly specified')
        measurement={'measurement_id':identity,'value':point['readout_percent'],'normalized_value':None,
            'assay_method':'MTT','unit':'percent non-treated control','representation':'summary',
            'asset':'pheno_data/Linn-2024/source/'+config['image'],'selector':f"Fig3b {point['cell']} {point['drug']} black point{point['point']} x={point['x_px']} y={point['y_px']}",
            'experimental_unit':None,'extraction':point,'reported_error_type':'SD','reported_technical_n':3,
            'independent_n':None,'readout_scale':'author-normalized summary, raw absorbance unavailable',
            'review':data['independent_review']}
        records.append({'comparison_id':identity,'source':source,'condition':cond,'phenotype':'cell_viability','assay_method':'MTT','batch':'Fig3b representative cytotoxicity experiment',
            'drug':point['drug'],'dose':None,'identity_evidence':'article.xml section4.1 names administered compound; exact salt/chemical form not reported; no cross-source identity merge attempted',
            'administered_form':'unspecified','dose_evidence':'Exact grid unpublished; pixel-calibrated log-axis estimate retained with15% uncertainty; canonical dose unresolved',
            'control_evidence':'Fig3 caption non-treated reference set100%; exact matched vehicle composition unavailable',
            'measurements':[measurement],'quality_reasons':['exact_source_dose_unresolved'],
            'statistics':{'treated_mean':point['readout_percent'],'control_mean':None,'scale_compatible':True,'control_matched':False,'p_raw':None,'p_adjusted':None,
                'evidence':'Caption defines non-treated normalization; matched vehicle not established. Technical-well SD cannot give independent-group p; CC50 independent n does not apply.'}})
    return records

def build_sources(selected):
    sources=[];records=[];audits={}
    for name in selected:
        context,data,config=source_context(name)
        sources.append(context)
        current=linn_records(name,data,config)
        records.extend(current)
        audits[name]={'curation_version':data['curation_version'],'source_comparisons':len(current),
            'independent_review':data['independent_review'],'status':'archived and fully extracted; exact dose grid unresolved; public build determines admission'}
    return {'schema_version':1,'sources':sources,'records':records,'source_audits':audits}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',choices=list(SOURCES),action='append')
    parser.add_argument('--prepare',action='store_true',help='explicitly freeze tracked reviewed curation into declared local additions')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    selected=args.source or list(SOURCES)
    if args.prepare:
        for name in selected: prepare(name)
    if args.output:
        args.output.write_text(json.dumps(build_sources(selected),indent=2,allow_nan=False)+'\n')
    elif not args.prepare: parser.error('--output or --prepare required')
