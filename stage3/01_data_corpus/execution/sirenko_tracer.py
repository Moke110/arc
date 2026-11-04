#!/usr/bin/env python3
"""Extract one original Sirenko comparison to demonstrate honest exclusion."""
import argparse
import hashlib
import json
from pathlib import Path

import xlrd


def extract(root):
    item='pheno_data/Sirenko-2013'
    asset=item+'/source/publisher-mmc3.xls'
    sheet=xlrd.open_workbook(root/asset).sheet_by_name('24hr_beat_count')
    if sheet.cell_value(1,3)!='Amoxicillin' or sheet.cell_value(1,2)!='26787-78-0':
        raise ValueError('Unexpected frozen Sirenko coordinate identity')
    assets=[]
    for line in (root/item/'manifest.sha256').read_text().splitlines():
        sha,name=line.split('  ',1)
        assets.append({'path':item+'/source/'+name,'sha256':sha})
    measurements=[]
    for replicate,raw_col,normalized_col in ((1,5,14),(2,6,15)):
        measurements.append({'measurement_id':f'sirenko_2013:24hr_beat_count:row2:dose0.3:r{replicate}',
            'value':sheet.cell_value(1,raw_col),'normalized_value':sheet.cell_value(1,normalized_col),
            'assay_method':'calcium_flux_peak_count','unit':'peaks in author acquisition window',
            'representation':'raw','asset':asset,'selector':f'24hr_beat_count!{chr(65+raw_col)}2; normalized {chr(65+normalized_col)}2',
            'experimental_unit':None,'replicate_kind':'technical_well_independence_unknown'})
    record={'comparison_id':'sirenko_2013:24hr_beat_count:row2:dose0.3','source':'sirenko_2013',
        'condition':{'condition_cell_type':'hiPSC-CM','condition_drug_time':24.,'condition_hypoxia_time':0.,
            'condition_reoxygen_time':0.,'condition_drug_schedule':'direct','condition_vehicle':'DMSO',
            'condition_species':'Homo sapiens','condition_source_protocol':'sirenko_2013_24h'},
        'phenotype':'beat_count','assay_method':'calcium_flux_peak_count','batch':None,'drug':None,'dose':.3e-6,
        'identity_evidence':None,'dose_evidence':'Supplemental Table 3 header concentration 0.3 micromolar',
        'control_evidence':None,'source_chemical_name':'Amoxicillin','source_cas':'26787-78-0',
        'measurements':measurements,'quality_reasons':['canonical_identity_and_pharmacology_not_yet_curated'],
        'statistics':{'scale_compatible':False,'control_matched':False,
            'evidence':'Source Control column is commercial classification, not matched vehicle; public vehicle values unavailable.'}}
    return {'schema_version':1,'sources':[{'source':'sirenko_2013','item':item,'test_only':False,
                'status':'archived','assets':assets}],'records':[record]}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    root=next(path for path in Path(__file__).resolve().parents if (path/'.git').exists())
    Path(args.output).write_text(json.dumps(extract(root),indent=2,sort_keys=True)+'\n')
