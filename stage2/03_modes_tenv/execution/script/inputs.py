#!/usr/bin/env python
"""Frozen tenv input release. Test truth is kept outside predictor manifests."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize


def repo_root():
    for path in Path(__file__).resolve().parents:
        if (path / '.git').exists():
            return path
    raise RuntimeError('repository root not found')

ROOT = repo_root()
STEP = Path(__file__).resolve().parent.parent.parent
DEFAULT_INPUT = STEP / 'output/input'
CASE_COLUMNS = ['case_id', 'drug_key', 'dose_um', 'control_index', 'group_id']
SEEDS = [17, 42, 2026]


def sha256(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')


def group_key(drug, dose):
    return hashlib.sha256(f'{drug}|{float(dose):.12g}'.encode()).hexdigest()[:24]


def normalized_structure(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f'invalid structure: {smiles}')
    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.FragmentParent(mol)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
    parent = Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)
    key = Chem.MolToInchiKey(mol)
    if not key:
        raise ValueError('parent InChIKey unavailable')
    return mol,parent,key.split('-')[0]


def split_drugs(drugs):
    """Shuffle structure blocks once, greedy balance three outer folds."""
    blocks = sorted(drugs.parent_connectivity.unique())
    rng = np.random.default_rng(20261105)
    rng.shuffle(blocks)
    size = drugs.groupby('parent_connectivity').size().to_dict()
    outer = [[],[],[]]
    counts = [0,0,0]
    for block in blocks:
        fold = int(np.argmin(counts))
        outer[fold].append(block)
        counts[fold] += size[block]
    records = []
    for fold in range(3):
        remaining = [b for b in blocks if b not in outer[fold]]
        val_rng = np.random.default_rng(20261105+fold+1)
        val_rng.shuffle(remaining)
        target = max(1,round(sum(size[b] for b in remaining)*.15))
        validation = []
        used = 0
        for block in remaining:
            if used < target:
                validation.append(block)
                used += size[block]
        for row in drugs.itertuples():
            role = 'test' if row.parent_connectivity in outer[fold] else ('validation' if row.parent_connectivity in validation else 'train')
            records.append(dict(fold=fold,drug_key=row.drug_key,parent_connectivity=row.parent_connectivity,role=role))
    return pd.DataFrame(records)


def _normalized_rows(matrix, rows, columns):
    # Dense h5py reads one short row slab; never construct observation × gene matrix.
    for start in range(0,len(rows),32):
        selected = rows[start:start+32]
        values = np.asarray(matrix[selected,:],dtype=np.float32)[:,columns]
        totals = values.sum(axis=1,keepdims=True)
        if np.any(totals <= 0):
            raise ValueError('zero projected library size')
        yield start,np.log1p(values*(np.float32(10000.)/totals)).astype(np.float32)


def _direction_labels(values,positive):
    eligible = np.flatnonzero(values > 0 if positive else values < 0)
    order = np.argsort(-values[eligible] if positive else values[eligible],kind='stable')[:100]
    labels = np.zeros(values.shape,dtype=bool)
    labels[eligible[order]] = True
    return labels


def build_release(corpus,genes_path,output,expected_genes=18164,expected_drugs=460):
    corpus,genes_path,output = map(Path,(corpus,genes_path,output))
    genes = pd.read_csv(genes_path,keep_default_na=False)
    if len(genes)!=expected_genes or genes.hgnc_id.duplicated().any():
        raise ValueError('wrong or duplicate gene axis')
    with h5py.File(corpus,'r') as source:
        obs,var = read_elem(source['obs']),read_elem(source['var'])
        obs = obs.reset_index(drop=True)
        var_lookup = {str(g):i for i,g in enumerate(var.hgnc_id)}
        if any(str(g) not in var_lookup for g in genes.hgnc_id):
            raise ValueError('missing gene identity')
        columns = np.array([var_lookup[str(g)] for g in genes.hgnc_id])
        blank = obs.is_blank.astype(bool).to_numpy()
        control = obs.is_control.astype(bool).to_numpy()
        dmso = control & ~blank & obs.compound_name.astype(str).str.upper().eq('DMSO').to_numpy() & obs.vehicle.astype(str).str.upper().eq('DMSO').to_numpy()
        treated = ~control & ~blank
        wells = obs.loc[treated].copy()
        if wells.well.duplicated().any() or wells.drug_key.isna().any() or not np.isfinite(wells.dose_um).all() or (wells.dose_um<=0).any():
            raise ValueError('invalid well/drug/dose identity')
        drugs = []
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
        fingerprints=[]
        for drug,rows in wells.groupby('drug_key',sort=True,observed=True):
            structures = [normalized_structure(s) for s in rows.smiles.astype(str).unique()]
            if len({s[1] for s in structures}) != 1:
                raise ValueError(f'inconsistent normalized parent: {drug}')
            mol,smiles,connectivity=structures[0]
            drugs.append(dict(drug_key=str(drug),parent_smiles=smiles,parent_connectivity=connectivity,morgan_index=len(drugs)))
            fingerprints.append(generator.GetFingerprintAsNumPy(mol))
        drugs=pd.DataFrame(drugs)
        if len(drugs)!=expected_drugs:
            raise ValueError(f'expected {expected_drugs} drugs, got {len(drugs)}')
        plates=sorted(wells.plate.astype(str).unique())
        controls=[]
        for plate in plates:
            rows=np.flatnonzero(dmso & obs.plate.astype(str).eq(plate).to_numpy())
            if not len(rows):
                raise ValueError(f'missing same-plate DMSO: {plate}')
            controls.append(dict(control_index=len(controls),plate=plate,n_wells=len(rows)))
        output.mkdir(parents=True,exist_ok=True)
        genes.to_csv(output/'genes.csv',index=False)
        drugs.to_csv(output/'drugs.csv',index=False)
        pd.DataFrame(controls).to_csv(output/'controls.csv',index=False)
        np.save(output/'morgan.npy',np.array(fingerprints,dtype=np.uint8))
        ctrl=np.lib.format.open_memmap(output/'controls.npy',mode='w+',dtype='float32',shape=(len(controls),len(genes)))
        for rec in controls:
            rows=np.flatnonzero(dmso & obs.plate.astype(str).eq(rec['plate']).to_numpy())
            total=np.zeros(len(genes),dtype=np.float64)
            for _,values in _normalized_rows(source['X'],rows,columns):
                total+=values.sum(axis=0,dtype=np.float64)
            ctrl[rec['control_index']]=total/len(rows)
        ctrl.flush()
        plate_indices={r['plate']:r['control_index'] for r in controls}
        cases=pd.DataFrame(dict(case_id=wells.well.astype(str).to_numpy(),drug_key=wells.drug_key.astype(str).to_numpy(),dose_um=wells.dose_um.to_numpy(),control_index=[plate_indices[str(p)] for p in wells.plate],group_id=[group_key(d,q) for d,q in zip(wells.drug_key,wells.dose_um)]))
        cases.to_csv(output/'wells.csv',index=False)
        folds=split_drugs(drugs)
        folds.to_csv(output/'folds.csv',index=False)
        source_rows=np.flatnonzero(treated)
        fold_counts=[]
        for fold in range(3):
            base=output/f'fold_{fold}'
            (base/'evaluator').mkdir(parents=True,exist_ok=True)
            membership=folds[folds.fold==fold]
            record={'fold':fold}
            for role in ['train','validation','test']:
                drug_set=set(membership[membership.role==role].drug_key)
                local_indices=np.flatnonzero(cases.drug_key.isin(drug_set))
                role_cases=cases.iloc[local_indices].copy().reset_index(drop=True)
                filename='test_cases.csv' if role=='test' else role+'.csv'
                role_cases.to_csv(base/filename,index=False)
                record[role]={'drugs':len(drug_set),'wells':len(role_cases),'groups':int(role_cases.group_id.nunique())}
                if role!='test':
                    destination=np.lib.format.open_memmap(base/f'{role}_expression.npy',mode='w+',dtype='float32',shape=(len(role_cases),len(genes)))
                    for start,values in _normalized_rows(source['X'],source_rows[local_indices],columns):
                        count=len(values)
                        destination[start:start+count]=values-ctrl[role_cases.control_index.to_numpy()[start:start+count]]
                    destination.flush()
                    del destination
                else:
                    group_order=list(dict.fromkeys(role_cases.group_id))
                    group_lookup={g:i for i,g in enumerate(group_order)}
                    sums=np.zeros((len(group_order),len(genes)),dtype=np.float64)
                    n=np.zeros(len(group_order),dtype=np.int64)
                    for start,values in _normalized_rows(source['X'],source_rows[local_indices],columns):
                        for offset,value in enumerate(values):
                            row=role_cases.iloc[start+offset]
                            j=group_lookup[row.group_id]
                            sums[j]+=value-ctrl[int(row.control_index)]
                            n[j]+=1
                    truth=(sums/n[:,None]).astype(np.float32)
                    np.save(base/'evaluator/truth.npy',truth)
                    np.save(base/'evaluator/positive.npy',np.array([_direction_labels(v,True) for v in truth]))
                    np.save(base/'evaluator/negative.npy',np.array([_direction_labels(v,False) for v in truth]))
                    groups=role_cases.drop_duplicates('group_id').set_index('group_id').loc[group_order,['drug_key','dose_um']].reset_index()
                    groups['n_wells']=n
                    groups.to_csv(base/'evaluator/groups.csv',index=False)
            write_json(base/'predictor.json',{'schema':'tenv-predictor-v1','fold':fold,'target_kind':'delta','gene_count':len(genes),'training_seeds':SEEDS,'files':{'genes':'../genes.csv','drugs':'../drugs.csv','morgan':'../morgan.npy','controls':'../controls.npy','control_metadata':'../controls.csv','train_cases':'train.csv','validation_cases':'validation.csv','test_cases':'test_cases.csv','train_expression':'train_expression.npy','validation_expression':'validation_expression.npy'}})
            fold_counts.append(record)
        del ctrl
    files={str(p.relative_to(output)):sha256(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    manifest={'schema':'tenv-input-v1','target_kind':'delta','gene_count':len(genes),'drug_count':len(drugs),'well_count':len(wells),'group_count':int(cases.group_id.nunique()),'split_seed':20261105,'training_seeds':SEEDS,'source_sha256':sha256(corpus),'gene_source_sha256':sha256(genes_path),'builder_sha256':sha256(Path(__file__)),'rdkit_version':rdBase.rdkitVersion,'morgan':{'radius':2,'bits':2048,'chirality':True,'binary':True,'invariants':'atom'},'normalization':'Cleanup; FragmentParent; Uncharger; canonical isomeric SMILES; connectivity=parent InChIKey first block','expression':'project raw counts onto ordered gene axis; CP10K; log1p; subtract same-plate DMSO mean log expression','row_counts':{'raw':len(obs),'treated':len(wells),'dmso':int(dmso.sum()),'other_controls':int((control&~blank&~dmso).sum()),'blank':int(blank.sum())},'folds':fold_counts,'files':files}
    write_json(output/'manifest.json',manifest)
    verify_release(output)
    return manifest


def predictor_manifest(path,fold):
    """Return only the frozen predictor contract, with safe absolute file paths."""
    path=Path(path)
    base=path/f'fold_{fold}'
    result=json.loads((base/'predictor.json').read_text())
    allowed={'genes','drugs','morgan','controls','control_metadata','train_cases','validation_cases','test_cases','train_expression','validation_expression'}
    if set(result['files'])!=allowed or result['target_kind']!='delta':
        raise ValueError('unsafe predictor manifest')
    expected={'genes':'../genes.csv','drugs':'../drugs.csv','morgan':'../morgan.npy','controls':'../controls.npy','control_metadata':'../controls.csv','train_cases':'train.csv','validation_cases':'validation.csv','test_cases':'test_cases.csv','train_expression':'train_expression.npy','validation_expression':'validation_expression.npy'}
    if result['files']!=expected:
        raise ValueError('unsafe predictor file path')
    result['files']={k:str((base/v).resolve()) for k,v in result['files'].items()}
    return result


def load_cases(path,fold,role):
    filename='test_cases.csv' if role=='test' else f'{role}.csv'
    if role not in {'train','validation','test'}:
        raise ValueError('unknown role')
    cases=pd.read_csv(Path(path)/f'fold_{fold}'/filename,keep_default_na=False)
    if list(cases.columns)!=CASE_COLUMNS or cases.case_id.duplicated().any():
        raise ValueError('unsafe or duplicate predictor cases')
    return cases


def verify_release(path):
    path=Path(path)
    manifest=json.loads((path/'manifest.json').read_text())
    for rel,digest in manifest['files'].items():
        file=path/rel
        if not file.resolve().is_relative_to(path.resolve()) or sha256(file)!=digest:
            raise ValueError(f'input hash mismatch: {rel}')
    drugs=pd.read_csv(path/'drugs.csv')
    folds=pd.read_csv(path/'folds.csv')
    if len(drugs)!=manifest['drug_count'] or len(pd.read_csv(path/'genes.csv'))!=manifest['gene_count']:
        raise ValueError('axis count mismatch')
    controls=pd.read_csv(path/'controls.csv')
    control_values=np.load(path/'controls.npy',mmap_mode='r')
    morgan=np.load(path/'morgan.npy',mmap_mode='r')
    if (control_values.shape!=(len(controls),manifest['gene_count'])
            or not np.isfinite(control_values).all()
            or controls.plate.duplicated().any()
            or not np.array_equal(controls.control_index,np.arange(len(controls)))
            or (controls.n_wells<=0).any()):
        raise ValueError('invalid control axis')
    if morgan.shape!=(len(drugs),2048) or not np.isin(morgan,[0,1]).all() or not np.array_equal(drugs.morgan_index,np.arange(len(drugs))):
        raise ValueError('invalid chemical axis')
    if not np.array_equal(np.unique(folds.fold),[0,1,2]):
        raise ValueError('missing folds')
    test_coverage=[]
    for fold in range(3):
        predictor_manifest(path,fold)
        roles={role:load_cases(path,fold,role) for role in ['train','validation','test']}
        membership=folds[folds.fold==fold]
        if len(membership)!=len(drugs) or membership.drug_key.duplicated().any() or set(membership.drug_key)!=set(drugs.drug_key):
            raise ValueError('incomplete drug membership')
        if membership.groupby('parent_connectivity').role.nunique().max()!=1:
            raise ValueError('parent-connectivity leakage')
        sets=[]
        for role,cases in roles.items():
            expected=set(membership[membership.role==role].drug_key)
            if set(cases.drug_key)!=expected:
                raise ValueError('case role leakage')
            if (not np.isfinite(cases.dose_um).all() or (cases.dose_um<=0).any()
                    or not cases.control_index.isin(controls.control_index).all()
                    or any(group_key(d,q)!=g for d,q,g in zip(cases.drug_key,cases.dose_um,cases.group_id))):
                raise ValueError('invalid case metadata')
            sets.append(set(cases.case_id))
            if role!='test':
                x=np.load(path/f'fold_{fold}/{role}_expression.npy',mmap_mode='r')
                if x.shape!=(len(cases),manifest['gene_count']) or any(not np.isfinite(x[i:i+32]).all() for i in range(0,len(x),32)):
                    raise ValueError('target shape mismatch')
        if sets[0]&sets[1] or sets[0]&sets[2] or sets[1]&sets[2]:
            raise ValueError('case leakage')
        test_coverage.extend(roles['test'].drug_key.unique())
        groups=pd.read_csv(path/f'fold_{fold}/evaluator/groups.csv')
        truth=np.load(path/f'fold_{fold}/evaluator/truth.npy',mmap_mode='r')
        if truth.shape!=(len(groups),manifest['gene_count']) or not np.isfinite(truth).all() or set(groups.group_id)!=set(roles['test'].group_id):
            raise ValueError('truth identity/shape mismatch')
        counts=roles['test'].groupby('group_id').size()
        if groups.group_id.duplicated().any() or any(int(counts[g])!=int(n) for g,n in zip(groups.group_id,groups.n_wells)):
            raise ValueError('group replicate mismatch')
        for positive,name in [(True,'positive'),(False,'negative')]:
            labels=np.load(path/f'fold_{fold}/evaluator/{name}.npy',mmap_mode='r')
            if labels.shape!=truth.shape or any(not np.array_equal(labels[i],_direction_labels(truth[i],positive)) for i in range(len(truth))):
                raise ValueError('directional truth label mismatch')
    if len(test_coverage)!=len(drugs) or set(test_coverage)!=set(drugs.drug_key):
        raise ValueError('outer test coverage mismatch')
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    build=sub.add_parser('build')
    build.add_argument('--corpus',type=Path,default=ROOT/'stage2/02_dataset_corpus/output/corpus/v2/tenv.h5ad')
    build.add_argument('--genes',type=Path,default=ROOT/'stage2/02_dataset_corpus/output/corpus/v2/tables/genes.csv')
    build.add_argument('--output',type=Path,default=DEFAULT_INPUT)
    build.add_argument('--expected-genes',type=int,default=18164)
    build.add_argument('--expected-drugs',type=int,default=460)
    verify=sub.add_parser('verify')
    verify.add_argument('--input',type=Path,default=DEFAULT_INPUT)
    args=parser.parse_args()
    manifest=build_release(args.corpus,args.genes,args.output,args.expected_genes,args.expected_drugs) if args.command=='build' else verify_release(args.input)
    print(json.dumps({'status':'verified','gene_count':manifest['gene_count'],'drug_count':manifest['drug_count'],'folds':manifest['folds']}))

if __name__=='__main__':
    main()
