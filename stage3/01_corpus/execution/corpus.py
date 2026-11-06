"""Build and verify source-first Heliste corpora using frozen local evidence only."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile

import openpyxl
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import rdBase

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
ITEM = ROOT/'pheno_data/Heliste-2021'
sys.path.insert(0, str(ITEM/'script'))
from desalt_chemicals import validate_record, POLICY

S1 = '41598_2021_96033_MOESM2_ESM.xlsx'
SOURCE_SHA = 'e3bce1a6b017d1d71a2033e89986b462a294fdf576bfeff16f9b75fe60c73aac'
HEADERS = ['Drug','CTRL_40nM','CTRL_200nM','CTRL_1000nM','CTRL_5000nM','HYP_40nM','HYP_200nM','HYP_1000nM','HYP_5000nM']
SCHEMA = pa.schema([
    ('source',pa.string()),('source_label',pa.string()),('canonical_smiles',pa.string()),
    ('dose_um',pa.float64()),('condition',pa.string()),('phenotype',pa.string()),
    ('source_asset',pa.string()),('source_cell',pa.string()),('rep',pa.int32()),
    ('readout',pa.struct([('source_reported',pa.float64())])),
    ('direction',pa.int8()),
])

DIRECTION_POLICY = {
    'version': 1,
    'method': 'magnitude_only',
    'reference': 'same-arm 0.01% DMSO vehicle',
    'reference_value': 100,
    'threshold_percent': 20,
    'zero_interval': [80,120],
    'zero_interval_inclusive': True,
    'normalization_status': 'working_assumption',
    'normalization_evidence': 'Paper confirms DMSO normalization; same-arm reference at 100 adopted by user, not explicitly confirmed in paper',
    'labels': {'-1': 'below lower boundary', '0': 'within magnitude band, not demonstrated equivalence', '1': 'above upper boundary'},
    'statistical_significance_claim': False,
}

def direction_from_readout(value):
    """Classify source relative viability under the frozen magnitude policy."""
    lower, upper = DIRECTION_POLICY['zero_interval']
    return -1 if value < lower else 1 if value > upper else 0

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def json_bytes(value):
    return (json.dumps(value,sort_keys=True,ensure_ascii=False,allow_nan=False,indent=2)+'\n').encode()

def frozen_inputs():
    """Validate archived bytes and individually declared local identity evidence."""
    log = json.loads((ITEM/'author_log.json').read_text())
    registry = json.loads((ROOT/'pheno_data/registry.json').read_text())
    if not any(r['item']=='Heliste-2021' and r['dataset_id']=='heliste_2021' and r['status']=='active' for r in registry['items']):
        raise ValueError('Heliste must be active in phenotype registry')
    inputs = {'pheno_data/registry.json':digest(ROOT/'pheno_data/registry.json')}
    declared = {entry['path']:entry['sha256'] for entry in log['sources']}
    source_manifest = dict((name,sha) for sha,name in (line.split(None,1) for line in (ITEM/'manifest.sha256').read_text().splitlines()))
    if declared != source_manifest or {p.name for p in (ITEM/'source').iterdir()} != set(declared):
        raise ValueError('source manifest/archive inventory mismatch')
    for entry in log['sources']:
        path = ITEM/'source'/entry['path']
        if digest(path) != entry['sha256'] or path.stat().st_size != entry['bytes']:
            raise ValueError('source provenance mismatch: '+str(path))
        inputs[str(path.relative_to(ROOT))] = digest(path)
    if digest(ITEM/'source'/S1) != SOURCE_SHA:
        raise ValueError('unsupported Dataset S1 bytes')
    additions = {e['path']:e for e in log['local_additions']}
    existing_local = {str(p.relative_to(ITEM)) for p in (ITEM/'local').rglob('*') if p.is_file()}
    if len(additions)!=len(log['local_additions']) or set(additions)!=existing_local:
        raise ValueError('undeclared or duplicate local evidence')
    identity_path = 'local/chemical_resolution_v2.json'
    if identity_path not in additions:
        raise ValueError('identity evidence must be declared in author_log')
    for name, entry in additions.items():
        path = ITEM/name
        if not name.startswith('local/') or not path.resolve().is_relative_to((ITEM/'local').resolve()):
            raise ValueError('local evidence outside item')
        if digest(path) != entry['sha256'] or path.stat().st_size != entry['bytes']:
            raise ValueError('local provenance mismatch: '+name)
        inputs[str(path.relative_to(ROOT))] = digest(path)
    for path in [ITEM/'author_log.json',ITEM/'manifest.sha256']:
        inputs[str(path.relative_to(ROOT))] = digest(path)
    resolution = json.loads((ITEM/identity_path).read_text())
    if resolution['original_snapshot_sha256'] != digest(ITEM/'local/chemical_resolution_v1.json'):
        raise ValueError('Original chemical reference snapshot mismatch')
    if resolution['canonicalization']['options'] != POLICY:
        raise ValueError('Unsupported desalted parent policy')
    if resolution['canonicalization']['version'] != rdBase.rdkitVersion:
        raise ValueError('RDKit version differs from frozen canonicalization')
    return resolution, inputs

def materialize(output):
    resolution, inputs = frozen_inputs()
    workbook = openpyxl.load_workbook(ITEM/'source'/S1,read_only=True,data_only=True)
    sheet = workbook.active
    rows = list(sheet.values)
    title = sheet.title
    workbook.close()
    if list(rows[3]) != HEADERS or len(rows[4:]) != 689:
        raise ValueError('unexpected Dataset S1 layout')
    labels = [row[0] for row in rows[4:]]
    identities = resolution['records']
    by_label = {r['source_label']:r for r in identities}
    if len(by_label)!=len(identities) or set(by_label)!=set(labels):
        raise ValueError('identity evidence does not exactly cover source labels')
    originals = {r['source_label']:r for r in json.loads((ITEM/'local/chemical_resolution_v1.json').read_text())['records']}
    for label, identity in by_label.items():
        validate_record(identity)
        if identity['original_reference'] != originals[label]:
            raise ValueError('Parent evidence differs from original frozen reference')
        if identity['status'] not in {'resolved','unresolved','ambiguous'}:
            raise ValueError('unknown resolution status')
        if identity['status']!='resolved' and not identity.get('reason'):
            raise ValueError('pending identity needs explicit reason')
    source_asset = str((ITEM/'source'/S1).relative_to(ROOT))
    output.mkdir(parents=True)
    coverage = {}
    for condition, start in [('baseline',1),('H-R',5)]:
        records, pending = [], []
        for excel_row,row in enumerate(rows[4:],5):
            identity = by_label[row[0]]
            cells = []
            for offset,dose in enumerate([0.04,0.2,1.0,5.0]):
                column = start+offset
                value = row[column]
                if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value):
                    raise ValueError('nonfinite source cell')
                cell = f'{title}!{openpyxl.utils.get_column_letter(column+1)}{excel_row}'
                cells.append(cell)
                if identity['status']=='resolved':
                    records.append(dict(source='heliste_2021',source_label=row[0],canonical_smiles=identity['canonical_smiles'],dose_um=dose,condition=condition,phenotype='cell_viability',source_asset=source_asset,source_cell=cell,rep=None,readout={'source_reported':float(value)},direction=direction_from_readout(value)))
            if identity['status']!='resolved':
                pending.append(dict(source_label=row[0],status=identity['status'],reason=identity['reason'],source_cells=cells,values=[float(row[start+i]) for i in range(4)]))
        directory = output/condition
        directory.mkdir()
        pq.write_table(pa.Table.from_pylist(records,schema=SCHEMA),directory/'corpus.parquet')
        with (directory/'corpus.csv').open('w',newline='') as stream:
            writer = csv.DictWriter(stream,fieldnames=SCHEMA.names,lineterminator="\n")
            writer.writeheader()
            for record in records:
                writer.writerow({**record,'rep':'','readout':json.dumps(record['readout'],sort_keys=True,allow_nan=False)})
        report = dict(source_label_rows=689,distinct_source_labels=len(set(labels)),resolved_labels=sum(r['status']=='resolved' for r in identities),unique_canonical_structures=len({r['canonical_smiles'] for r in identities if r['status']=='resolved'}),source_values=2756,admitted_values=len(records),pending_values=len(pending)*4,pending=pending,chemical_resolution=identities,canonicalization=resolution['canonicalization'])
        annotation = dict(source='Heliste 2021 Dataset S1',doi='10.1038/s41598-021-96033-z',source_asset=source_asset,condition=condition,source_arm='CTRL' if condition=='baseline' else 'HYP',phenotype='cell_viability',assay='CellTiter-Glo (CTG), ATP-dependent luminescence',cell_model='HL-1 mouse immortalized atrial cardiomyocytes',pretreatment_hours=24,hypoxia_hours=0 if condition=='baseline' else 48,reoxygenation_hours=0 if condition=='baseline' else 24,elapsed_hours_approximately=96,normoxia='duration-matched parallel normoxia' if condition=='baseline' else None,vehicle='0.01% DMSO',scale='author-supplied normalized relative viability; preserved exactly',normalization_denominator='exact oxygen-specific denominator unresolved',drug_persistence='unresolved',replicate_identity=None,measurement_semantics='source_reported; individual versus group mean unknown',absent_evidence=['raw RLU','vehicle measurements','replicate identities','SD','n','p-values','source endpoint labels'],admission='resolved chemical identity; no endpoint or significance gate',dose_unit='micromolar',direction_policy=DIRECTION_POLICY,canonicalization=POLICY,evidence=['PMC8371130.1.xml','41598_2021_96033_MOESM1_ESM.pdf',S1])
        (directory/'chemical_resolution_report.json').write_bytes(json_bytes(report))
        (directory/'source_annotation.json').write_bytes(json_bytes(annotation))
        coverage[condition] = {k:report[k] for k in ['source_label_rows','resolved_labels','unique_canonical_structures','source_values','admitted_values','pending_values']}
    manifest = dict(schema_version=2,inputs=inputs,coverage=coverage,code={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__).resolve(),ITEM/'script/resolve_chemicals.py',ITEM/'script/desalt_chemicals.py']},software={'pyarrow':pa.__version__,'rdkit':rdBase.rdkitVersion,'openpyxl':openpyxl.__version__},files={str(p.relative_to(output)):digest(p) for p in sorted(output.rglob('*')) if p.is_file()})
    (output/'manifest.json').write_bytes(json_bytes(manifest))
    return manifest

def verify(output):
    manifest = json.loads((output/'manifest.json').read_text())
    expected_files = set(manifest['files'])|{'manifest.json'}
    if {str(p.relative_to(output)) for p in output.rglob('*') if p.is_file()}!=expected_files:
        raise ValueError('unexpected or missing release files')
    for name, sha in manifest['files'].items():
        path = output/name
        if not path.resolve().is_relative_to(output.resolve()) or digest(path)!=sha:
            raise ValueError('output checksum mismatch: '+name)
    with tempfile.TemporaryDirectory() as temp:
        rebuilt = Path(temp)/'rebuilt'
        current = materialize(rebuilt)
        if current!=manifest:
            raise ValueError('release differs from offline rebuild or current provenance')
        for condition in ['baseline','H-R']:
            actual = pq.read_table(output/condition/'corpus.parquet')
            if actual.schema!=SCHEMA or not actual.equals(pq.read_table(rebuilt/condition/'corpus.parquet')):
                raise ValueError('Parquet source fidelity failure')
            if (output/condition/'corpus.csv').read_bytes()!=(rebuilt/condition/'corpus.csv').read_bytes():
                raise ValueError('CSV source fidelity failure')
    return manifest

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['build','verify'])
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent.parent/'output')
    args = parser.parse_args()
    output = args.output.resolve()
    if args.command=='verify' or output.exists():
        manifest = verify(output)
    else:
        output.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=output.parent) as temp:
            staging = Path(temp)/'release'
            manifest = materialize(staging)
            staging.rename(output)
    print(json.dumps(manifest['coverage'],sort_keys=True))

if __name__=='__main__':
    try:
        main()
    except (ValueError,KeyError,OSError) as error:
        print(str(error),file=sys.stderr)
        sys.exit(1)
