"""Heliste Dataset S1: preserve relative CTG values without inventing controls."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import math
import xml.etree.ElementTree as ET
import zipfile
import openpyxl

S1 = '41598_2021_96033_MOESM2_ESM.xlsx'
S3 = '41598_2021_96033_MOESM4_ESM.xlsx'


def _independent_counts(path):
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as archive:
        doc = ET.fromstring(archive.read('xl/worksheets/sheet1.xml'))
        cells = [c for row in doc.findall('s:sheetData/s:row', ns) if int(row.attrib['r']) >= 5 for c in row.findall('s:c', ns)]
        return Counter(c.attrib['r'].rstrip('0123456789') for c in cells if c.attrib.get('t') not in {'s','inlineStr'} and c.find('s:v', ns) is not None)


def build_bundle(root):
    root = Path(root)
    item = root/'pheno_data/Heliste-2021'
    source_dir = item/'source'
    rows = list(openpyxl.load_workbook(source_dir/S1, read_only=True, data_only=True).active.values)
    expected_headers = ['Drug','CTRL_40nM','CTRL_200nM','CTRL_1000nM','CTRL_5000nM','HYP_40nM','HYP_200nM','HYP_1000nM','HYP_5000nM']
    if list(rows[3]) != expected_headers:
        raise ValueError('Heliste S1 header changed; source adapter review required')
    data = rows[4:]
    labels = [r[0] for r in data]
    if any(not isinstance(v,str) or not v.strip() for v in labels):
        raise ValueError('invalid S1 drug label')
    counts = _independent_counts(source_dir/S1)
    finite = [sum(isinstance(r[i],(float,int)) and not isinstance(r[i],bool) and math.isfinite(r[i]) for r in data) for i in range(1,9)]
    if any(counts[letter] != count for letter,count in zip('BCDEFGHI',finite)):
        raise ValueError('independent OOXML/openpyxl numeric counts disagree')
    curation = item/'local/assessment_v1.json'
    frozen_assessment = json.loads(curation.read_text())
    source_assets = [{'path':str(path.relative_to(root)), 'sha256':hashlib.sha256(path.read_bytes()).hexdigest()} for path in sorted([*source_dir.iterdir(), curation])]
    records = []
    for excel_row, row in enumerate(data,5):
        for column, header in enumerate(expected_headers[1:],1):
            oxygen, dose = header.split('_')
            is_hr = oxygen == 'HYP'
            value = row[column]
            if not isinstance(value,(float,int)) or not math.isfinite(value):
                raise ValueError('nonfinite Heliste value: '+str(excel_row)+header)
            condition = {
                'condition_cell_type':'HL-1 mouse immortalized atrial cardiomyocyte',
                'condition_species':'Mus musculus',
                'condition_drug_time':24.0 if is_hr else None,
                'condition_hypoxia_time':48.0 if is_hr else 0.0,
                'condition_reoxygen_time':24.0 if is_hr else 0.0,
                'condition_drug_schedule':'24h pretreatment then HR; drug persistence unknown' if is_hr else 'parallel normoxia; drug persistence unknown',
                'condition_elapsed_time':96.0,
                'condition_drug_persistence':'unknown',
                'condition_vehicle':'0.01% DMSO',
                'condition_oxygen_protocol':'1% O2 hypoxia then reoxygenation' if is_hr else '21% O2 normoxia',
                'condition_source_context':'heliste_2021 drug library screen',
            }
            comparison = f'heliste_2021:S1:{excel_row}:{header}'
            records.append({
                'comparison_id':comparison,'source':'heliste_2021','condition':condition,
                'phenotype':'ctg_atp_viability','assay_method':'CTG','batch':'unresolved:heliste_2021:S1',
                'drug':None,'administered_drug_label':row[0], 'dose':int(dose[:-2])*1e-9,
                'dose_original':int(dose[:-2]),'dose_original_unit':'nM',
                'identity_evidence':'Dataset S1 drug label only; no chemical identifiers; aliases, salts and active parent unresolved',
                'dose_evidence':f'Dataset S1 column {header}; canonical mol/L = nM * 1e-9',
                'control_evidence':'Methods normalize to DMSO-treated controls; oxygen-specific denominator and actual control values not supplied',
                'quality_reasons':['unresolved_active_parent_identity','unresolved_ctg_normalization_denominator','independent_experimental_unit_unknown'] + ([] if is_hr else ['parallel_normoxia_elapsed96h_outside_direct24_48h_target']),
                'measurements':[{'measurement_id':comparison+':value','value':float(value),'assay_method':'CTG','unit':'relative_viability_supplied_scale','representation':'normalized_only','asset':str((source_dir/S1).relative_to(root)),'selector':f'Sheet1!{chr(65+column)}{excel_row}','experimental_unit':None,'normalized_value':None}],
                'statistics':{'control_matched':False,'scale_compatible':False,'family_complete':False,'family_matches':False,'evidence':'No DMSO values, normalization denominator, raw RLU, replicate identifiers, n, error or numerical p provided in S1'},
            })
    hitrows = list(openpyxl.load_workbook(source_dir/S3, read_only=True, data_only=True).active.values)
    audit = {'dataset_id':'heliste_2021','table':'Supplementary Dataset S1','source_labels':len(labels),'distinct_exact_labels':len(set(labels)),'numeric_values_by_oxygen':{'CTRL':sum(finite[:4]),'HYP':sum(finite[4:])},'numeric_values_by_column':dict(zip(expected_headers[1:],finite)),'independent_ooxml_column_counts':dict(counts),'dataset_s2':'shRNA; excluded from drug measurements','dataset_s3_hit_labels':sum(bool(r[0]) for r in hitrows[3:]),'canonical_drugs_resolved':0,'source_units':'relative viability, not raw RLU','normalization_denominator':'unresolved','corpus_expected_status':'empty-after-QC','view500_assessment':'not run','model_training':'not run'}
    for key in ('distinct_exact_labels','numeric_values_by_oxygen','dataset_s3_hit_labels'):
        if audit[key] != frozen_assessment[key]:
            raise ValueError('Frozen Heliste assessment differs: '+key)
    audit['curation_asset'] = str(curation.relative_to(root))
    audit['methods_evidence'] = frozen_assessment['methods_evidence']
    audit['ancillary_not_archived'] = frozen_assessment['ancillary_not_archived']
    return {'schema_version':1,'sources':[{'source':'heliste_2021','item':'pheno_data/Heliste-2021','status':'active','test_only':False,'assets':source_assets}],'records':records}, audit


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root = next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
    bundle, audit = build_bundle(root)
    bundle['source_audits'] = {'heliste_2021': audit}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(bundle, sort_keys=True, ensure_ascii=False, allow_nan=False, indent=2)+'\n')
    print(json.dumps(audit,sort_keys=True))

if __name__ == '__main__':
    main()
